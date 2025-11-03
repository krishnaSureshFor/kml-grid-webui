import os
import zipfile
import tempfile
from flask import Flask, render_template, request, send_file, redirect, url_for, flash
from fastkml import kml
from shapely.geometry import shape, Polygon, MultiPolygon, box, Point, GeometryCollection
from shapely.ops import transform, unary_union
import simplekml
from pyproj import Transformer, CRS
import zipfile

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
    Parses KML or KMZ file specified by filepath; returns list of shapely Polygon objects
    extracted from all Placemarks in the file, including nested MultiGeometry.
    """
    geometries = []

    # Read .kml or .kmz
    ext = filepath.lower().rsplit(".", 1)[-1]
    if ext == "kmz":
        with zipfile.ZipFile(filepath, "r") as z:
            kml_files = [f for f in z.namelist() if f.lower().endswith(".kml")]
            if not kml_files:
                return geometries
            with z.open(kml_files[0]) as f:
                doc = f.read()
    else:
        with open(filepath, "rb") as f:
            doc = f.read()

    k = kml.KML()
    k.from_string(doc)

    def extract_feats(feats):
        for feat in feats:
            # If the feature has geometry
            if hasattr(feat, "geometry") and feat.geometry:
                try:
                    geom = shape(feat.geometry)
                except Exception:
                    geom = None
                if geom:
                    # If geometry is a polygon
                    if isinstance(geom, Polygon):
                        geometries.append(geom)
                    # If MultiPolygon or collection, flatten it
                    elif isinstance(geom, MultiPolygon):
                        for g in geom.geoms:
                            if isinstance(g, Polygon):
                                geometries.append(g)
                    elif hasattr(geom, 'geoms'):
                        # For GeometryCollection etc
                        for sub in geom.geoms:
                            if isinstance(sub, Polygon):
                                geometries.append(sub)

            # Recurse sub-features
            if hasattr(feat, "features"):
                sub = feat.features
                if isinstance(sub, list) and sub:
                    extract_feats(sub)

    top_feats = k.features
    if isinstance(top_feats, list):
        extract_feats(top_feats)

    # Optional: unify geometries if you want one combined polygon
    # unified = unary_union(geometries)
    # if isinstance(unified, Polygon):
    #    return [unified]

    return geometries

def make_grid_for_polygon(poly, grid_meters, id_start, prefix="G"):
    """
    Create grid cells (clipped to polygon) in metric CRS and return geometries back to lon/lat.
    Returns list of tuples: (clipped_polygon_lonlat, center_point_lonlat, id_str), and next index.
    """
    # project to Web Mercator (meters)
    crs_from = CRS.from_epsg(4326)
    crs_to = CRS.from_epsg(3857)
    project_to_m = Transformer.from_crs(crs_from, crs_to, always_xy=True).transform
    project_to_lonlat = Transformer.from_crs(crs_to, crs_from, always_xy=True).transform

    poly_m = transform(project_to_m, poly)
    minx, miny, maxx, maxy = poly_m.bounds

    cells = []
    # compute counts (ensure at least 1)
    nx = max(1, int((maxx - minx) / grid_meters) + 1)
    ny = max(1, int((maxy - miny) / grid_meters) + 1)
    idx = id_start

    for i in range(nx):
        for j in range(ny):
            x0 = minx + i * grid_meters
            y0 = miny + j * grid_meters
            cell = box(x0, y0, x0 + grid_meters, y0 + grid_meters)
            if poly_m.intersects(cell):
                inter = poly_m.intersection(cell)
                # transform clipped geometry back to lon/lat
                cell_lonlat = transform(project_to_lonlat, inter)
                center_m = cell.centroid
                center_lonlat = transform(project_to_lonlat, center_m)
                id_str = f"{prefix}_{idx:04d}"
                cells.append((cell_lonlat, Point(center_lonlat.x, center_lonlat.y), id_str))
                idx += 1

    return cells, idx

@app.route("/", methods=["GET", "POST"])
def index():
    print("\n---- NEW REQUEST ----")

    if request.method == "POST":
        print("POST received")

        # Read form
        file = request.files.get("kmlfile")
        grid_size_raw = request.form.get("grid_size")
        prefix = request.form.get("prefix") or "G"
        preview = request.form.get("preview") == "on"

        print("File object:", file)
        print("Filename:", file.filename if file else None)
        print("Grid size raw:", grid_size_raw)
        print("Prefix:", prefix)
        print("Preview flag:", preview)

        # If file missing
        if not file or file.filename == "":
            print("ERROR: No file received")
            flash("File upload FAILED — server did NOT receive any file.", "danger")
            return redirect(request.url)

        # Grid size convert
        try:
            grid_size = float(grid_size_raw)
            print("Grid size OK:", grid_size)
        except Exception as e:
            print("ERROR: grid_size conversion failed:", e)
            flash("Grid size invalid: " + str(e), "danger")
            return redirect(request.url)

        # Save uploaded file
        try:
            tmpdir = tempfile.mkdtemp(dir=app.config["UPLOAD_FOLDER"])
            filepath = os.path.join(tmpdir, file.filename)
            file.save(filepath)
            print("File saved to:", filepath)
        except Exception as e:
            print("ERROR saving file:", e)
            flash("Server failed to save uploaded file: " + str(e), "danger")
            return redirect(request.url)

        # Parse KML polygons
        try:
            polys = parse_kml(filepath)
            print("Polygons found:", len(polys))
        except Exception as e:
            print("ERROR parsing KML:", e)
            flash("KML parsing failed: " + str(e), "danger")
            return redirect(request.url)

        if not polys:
            print("ERROR: No polygon geometries detected")
            flash("Uploaded KML contains NO polygons!", "danger")
            return redirect(request.url)

        # Generate grid
        try:
            print("Starting grid generation...")
            kml_out = simplekml.Kml()
            cur_idx = 1
            total_cells = 0

            for i, poly in enumerate(polys, start=1):
                cells, cur_idx = make_grid_for_polygon(poly, grid_size, cur_idx, prefix)
                total_cells += len(cells)
                print(f"Polygon {i}: cells = {len(cells)}")

                folder = kml_out.newfolder(name=f"Polygon_{i}")
                for poly_geom, center_point, nid in cells:
                    pg = folder.newpolygon(
                        name=nid,
                        outerboundaryis=list(poly_geom.exterior.coords)
                    )
                    pg.style.polystyle.fill = 0
                    pg.style.linestyle.color = simplekml.Color.red

                    pt = folder.newpoint(name=nid, coords=[(center_point.x, center_point.y)])
                    pt.style.iconstyle.scale = 0

            out_path = os.path.join(tmpdir, "Grid.kml")
            kml_out.save(out_path)
            print("Grid saved to:", out_path)

        except Exception as e:
            print("ERROR generating grid:", e)
            flash("Grid generation failed: " + str(e), "danger")
            return redirect(request.url)

        # Preview mode
        if preview:
            preview_url = url_for("download_file")
            print("Preview URL:", preview_url)
            return render_template(
                "index.html",
                preview_url=preview_url,
                download_url=preview_url,
                total_cells=total_cells,
                grid_size=grid_size,
                prefix=prefix
            )

        return send_file(out_path, as_attachment=True)

    # GET request → show form
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

