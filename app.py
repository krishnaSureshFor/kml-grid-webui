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
    Full recursive KML/KMZ parser.
    Supports: Document, Folder, Placemark, MultiGeometry, Polygon.
    Returns: list of shapely polygons in WGS84.
    """

    import zipfile
    from fastkml import kml
    from shapely.geometry import shape, Polygon, MultiPolygon

    geometries = []

    # --- READ KML/KMZ FILE ---
    ext = filepath.lower().rsplit(".", 1)[1]
    if ext == "kmz":
        with zipfile.ZipFile(filepath, "r") as z:
            kml_files = [n for n in z.namelist() if n.lower().endswith(".kml")]
            if not kml_files:
                return geometries
            with z.open(kml_files[0]) as f:
                doc = f.read()
    else:
        with open(filepath, "rb") as f:
            doc = f.read()

    k = kml.KML()
    k.from_string(doc)

    # --- FIX: features is a list, not callable ---
    def walk(feat_list):
        for feat in feat_list:
            # --- If it has geometry ---
            if hasattr(feat, "geometry") and feat.geometry:
                try:
                    g = shape(feat.geometry)

                    # Polygon
                    if isinstance(g, Polygon):
                        geometries.append(g)

                    # MultiPolygon
                    elif isinstance(g, MultiPolygon):
                        for sub in g.geoms:
                            geometries.append(sub)

                except Exception:
                    pass

            # --- Recurse into sub-features ---
            if hasattr(feat, "features"):
                sub = feat.features   # not callable
                if isinstance(sub, list):
                    walk(sub)

    top = k.features
    if isinstance(top, list):
        walk(top)

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
          @app.route("/", methods=("GET", "POST"))
@app.route("/", methods=("GET", "POST"))
def index():
    if request.method == "POST":

        # ---- READ FORM ----
        file = request.files.get("kmlfile")
        grid_size = request.form.get("grid_size")
        prefix = request.form.get("prefix") or "G"
        preview = request.form.get("preview") == "on"

        # ---- FILE CHECK ----
        if not file or file.filename == "":
            flash("Please upload a KML/KMZ file.", "danger")
            return redirect(request.url)

        if not allowed_file(file.filename):
            flash("Only .kml or .kmz allowed.", "danger")
            return redirect(request.url)

        # ---- GRID SIZE ----
        try:
            grid_size = float(grid_size)
        except:
            flash("Grid size must be a number.", "danger")
            return redirect(request.url)

        # ---- SAVE UPLOADED FILE ----
        tmpdir = tempfile.mkdtemp(dir=app.config["UPLOAD_FOLDER"])
        filepath = os.path.join(tmpdir, file.filename)
        file.save(filepath)

        # ---- PARSE KML ----
        polys = parse_kml(filepath)
        if not polys:
            flash("No polygon found in KML file.", "danger")
            return redirect(request.url)

        # ---- GENERATE GRID ----
        kml_out = simplekml.Kml()
        cur_idx = 1
        total_cells = 0

        for i, poly in enumerate(polys, start=1):
            cells, cur_idx = make_grid_for_polygon(poly, grid_size, cur_idx, prefix)
            total_cells += len(cells)

            folder = kml_out.newfolder(name=f"Polygon_{i}")
            for poly_geom, center_point, nid in cells:
                pg = folder.newpolygon(
                    name=nid,
                    outerboundaryis=list(poly_geom.exterior.coords)
                )
                pg.style.polystyle.fill = 0
                pg.style.linestyle.color = simplekml.Color.red

                # label point
                pt = folder.newpoint(name=nid, coords=[(center_point.x, center_point.y)])
                pt.style.iconstyle.scale = 0  # hide icon

        out_path = os.path.join(tmpdir, "Grid.kml")
        kml_out.save(out_path)

        # ---- IF PREVIEW ----
        if preview:
            # This creates a URL like: /download/somefolder/Grid.kml
            preview_url = url_for("download_file", path=os.path.basename(tmpdir) + "/Grid.kml")

            return render_template(
                "index.html",
                preview_url=preview_url,
                download_url=preview_url,
                total_cells=total_cells,
                grid_size=grid_size,
                prefix=prefix
            )

        # ---- ELSE DIRECT DOWNLOAD ----
        return send_file(out_path, as_attachment=True)

    # ---- INITIAL FORM ----
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
