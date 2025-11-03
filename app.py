import os
import io
import zipfile
import tempfile
from pathlib import Path
from flask import Flask, render_template, request, send_file, redirect, url_for, flash
from fastkml import kml
from shapely.geometry import shape, box, mapping, Point
from shapely.ops import transform
import simplekml
from pyproj import Transformer, CRS

# Configuration
UPLOAD_FOLDER = "uploads"
ALLOWED_EXT = {"kml", "kmz"}
os.makedirs(UPLOAD_FOLDER, exist_ok=True)

app = Flask(__name__)
app.secret_key = "change_me_to_a_secure_random_value"
app.config["UPLOAD_FOLDER"] = UPLOAD_FOLDER
app.config["MAX_CONTENT_LENGTH"] = 20 * 1024 * 1024  # 20 MB limit

def allowed_file(filename):
    return "." in filename and filename.rsplit(".", 1)[1].lower() in ALLOWED_EXT

def parse_kml(filepath):
    """
    Parse a KML (or KMZ) and return list of shapely polygons (in WGS84 lon/lat).
    Accepts KML with multiple Placemarks; extracts Polygon & MultiPolygon geometries.
    """
    geometries = []
    ext = filepath.lower().rsplit(".", 1)[1]
    if ext == "kmz":
        # Fast approach: unzip and find .kml inside
        with zipfile.ZipFile(filepath, 'r') as z:
            # find first .kml
            kml_names = [n for n in z.namelist() if n.lower().endswith(".kml")]
            if not kml_names:
                return geometries
            with z.open(kml_names[0]) as kf:
                doc = kf.read()
    else:
        with open(filepath, "rb") as f:
            doc = f.read()
    k = kml.KML()
    k.from_string(doc)
    # traverse features recursively
    def recurse_features(feats):
        for feat in feats:
            if hasattr(feat, "geometry") and feat.geometry:
                geom = feat.geometry
                try:
                    s = shape(geom)
                    geometries.append(s)
                except Exception:
                    pass
            if hasattr(feat, "features"):
                recurse_features(list(feat.features()))
    recurse_features(list(k.features()))
    return geometries

def make_grid_for_polygon(poly, grid_meters, id_start, prefix="G"):
    """
    poly: shapely polygon in WGS84 lon/lat
    grid_meters: cell size in meters
    id_start: integer starting index for numbering
    returns list of tuples: (cell_polygon_in_lonlat, center_point_in_lonlat, id_str)
    """
    # Reproject to an equal-area / metric projection for grid creation.
    # We'll use Web Mercator (EPSG:3857) for simplicity — accurate for moderate areas.
    crs_from = CRS.from_epsg(4326)
    crs_to = CRS.from_epsg(3857)
    project_to_m = Transformer.from_crs(crs_from, crs_to, always_xy=True).transform
    project_to_lonlat = Transformer.from_crs(crs_to, crs_from, always_xy=True).transform

    poly_m = transform(project_to_m, poly)
    minx, miny, maxx, maxy = poly_m.bounds

    # build grid in metric CRS
    cells = []
    nx = int((maxx - minx) / grid_meters) + 1
    ny = int((maxy - miny) / grid_meters) + 1
    idx = id_start
    for i in range(nx):
        for j in range(ny):
            x0 = minx + i * grid_meters
            y0 = miny + j * grid_meters
            cell = box(x0, y0, x0 + grid_meters, y0 + grid_meters)
            if poly_m.intersects(cell):
                inter = poly_m.intersection(cell)
                # Use the full cell polygon (not necessarily intersection) if any intersection wanted.
                # We'll use the clipped polygon (intersection) for cleaner KML coverage.
                cell_lonlat = transform(project_to_lonlat, inter)
                center_m = cell.centroid
                center_lonlat = transform(project_to_lonlat, center_m)
                id_str = f"{prefix}_{idx:04d}"
                cells.append((cell_lonlat, Point(center_lonlat.x, center_lonlat.y), id_str))
                idx += 1
    return cells, idx

@app.route("/", methods=("GET", "POST"))
def index():
    if request.method == "POST":
        # file upload
        file = request.files.get("kmlfile")
        grid_size = request.form.get("grid_size")
        prefix = request.form.get("prefix") or "G"
        preview = request.form.get("preview") == "on"
        if not file or file.filename == "":
            flash("Please upload a KML or KMZ file.", "danger")
            return redirect(request.url)
        if not allowed_file(file.filename):
            flash("Only .kml or .kmz files are allowed.", "danger")
            return redirect(request.url)
        try:
            grid_size = float(grid_size)
            assert grid_size > 0
        except Exception:
            flash("Please enter a valid grid size (in meters).", "danger")
            return redirect(request.url)

        # save uploaded file temporarily
        tmpdir = tempfile.mkdtemp(dir=app.config["UPLOAD_FOLDER"])
        filepath = os.path.join(tmpdir, file.filename)
        file.save(filepath)

        # parse KML to polygons
        polys = parse_kml(filepath)
        if not polys:
            flash("No polygon geometry found in uploaded KML/KMZ.", "danger")
            return redirect(request.url)

        # create KML
        kml_out = simplekml.Kml()
        # track numbering across polygons
        cur_idx = 1
        total_cells = 0

        for p_index, poly in enumerate(polys, start=1):
            cells, cur_idx = make_grid_for_polygon(poly, grid_size, cur_idx, prefix=prefix)
            total_cells += len(cells)
            # group in KML folder per polygon
            folder = kml_out.newfolder(name=f"Polygon_{p_index}_grids")
            for cell_poly, center_point, id_str in cells:
                pol = folder.newpolygon(name=id_str, outerboundaryis=list(cell_poly.exterior.coords))
                pol.style.polystyle.fill = 0  # no fill
                pol.style.linestyle.width = 1
                pol.style.linestyle.color = simplekml.Color.hex("ff0000")  # red
                # add point label as placemark (center)
                p = folder.newpoint(name=id_str, coords=[(center_point.x, center_point.y)])
                p.style.labelstyle.scale = 0.9
                p.style.iconstyle.scale = 0.0  # hide icon, keep label
            # optionally include original polygon outline as reference
            orig = folder.newpolygon(name=f"Polygon_{p_index}_source", outerboundaryis=list(poly.exterior.coords))
            orig.style.polystyle.fill = 0
            orig.style.linestyle.color = simplekml.Color.hex("ff00ff00")
        # add metadata placemark or description
        kml_out.document.description = f"Grid generated by web UI. Developed by Krishnasureshkumar FG"

        # save KML to a temporary file and send to user
        out_kml_path = os.path.join(tmpdir, "Grid.kml")
        kml_out.save(out_kml_path)

        # If preview requested, return page with embedded preview (leaflet + omnivore)
        if preview:
            # create a URL to serve the generated KML via Flask static send_file route
            # We'll provide a simple link to download and preview in the page.
            return render_template(
                "index.html",
                preview_url=url_for("download_file", path=os.path.basename(tmpdir) + "/" + "Grid.kml"),
                download_url=url_for("download_file", path=os.path.basename(tmpdir) + "/" + "Grid.kml"),
                total_cells=total_cells,
                grid_size=grid_size,
                prefix=prefix
            )
        else:
            return send_file(out_kml_path, as_attachment=True, download_name="Grid.kml")
    return render_template("index.html")

@app.route("/download/<path:path>")
def download_file(path):
    # path is relative path inside UPLOAD_FOLDER
    full = os.path.join(app.config["UPLOAD_FOLDER"], path)
    if os.path.exists(full):
        return send_file(full, as_attachment=True)
    else:
        flash("File not found.", "danger")
        return redirect(url_for("index"))

if __name__ == "__main__":
    app.run(debug=True, port=5000)
