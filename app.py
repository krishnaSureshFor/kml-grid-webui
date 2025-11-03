import os
import tempfile
from pathlib import Path
from flask import Flask, render_template, request, send_file, redirect, url_for, flash
from shapely.geometry import Polygon, Point, box
from shapely.ops import transform
import simplekml
from pyproj import Transformer, CRS
from lxml import etree

# ----------------------------
# Flask setup
# ----------------------------
UPLOAD_FOLDER = "uploads"
os.makedirs(UPLOAD_FOLDER, exist_ok=True)

app = Flask(__name__)
app.secret_key = "secret123"
app.config["UPLOAD_FOLDER"] = UPLOAD_FOLDER
app.config["MAX_CONTENT_LENGTH"] = 20 * 1024 * 1024  # 20 MB

ALLOWED_EXT = {"kml", "kmz"}


# ----------------------------
# Helper functions
# ----------------------------

def allowed_file(filename):
    return "." in filename and filename.rsplit(".", 1)[1].lower() in ALLOWED_EXT


# ✅ Manual XML KML parser (works for ANY polygon format)
def parse_kml(filepath):
    with open(filepath, "rb") as f:
        xml = f.read()

    root = etree.fromstring(xml)

    ns = {"kml": "http://www.opengis.net/kml/2.2"}
    polygons = []

    coords_list = root.findall(".//kml:Polygon//kml:coordinates", ns)

    for coords in coords_list:
        raw = coords.text.strip()

        # normalize
        parts = raw.replace("\n", " ").replace("\t", " ").split()

        ring = []
        for p in parts:
            vals = p.split(",")
            if len(vals) >= 2:
                lon = float(vals[0])
                lat = float(vals[1])
                ring.append((lon, lat))

        if len(ring) >= 3:
            polygons.append(Polygon(ring))

    return polygons


# ✅ Grid generator
def make_grid_for_polygon(poly, grid_m, id_start, prefix="G"):
    crs_from = CRS.from_epsg(4326)
    crs_to = CRS.from_epsg(3857)
    project_to_m = Transformer.from_crs(crs_from, crs_to, always_xy=True).transform
    project_to_ll = Transformer.from_crs(crs_to, crs_from, always_xy=True).transform

    poly_m = transform(project_to_m, poly)
    minx, miny, maxx, maxy = poly_m.bounds

    cells = []
    idx = id_start

    nx = int((maxx - minx) / grid_m) + 1
    ny = int((maxy - miny) / grid_m) + 1

    for i in range(nx):
        for j in range(ny):
            x0 = minx + i * grid_m
            y0 = miny + j * grid_m

            cell = box(x0, y0, x0 + grid_m, y0 + grid_m)

            if poly_m.intersects(cell):
                inter = poly_m.intersection(cell)
                inter_ll = transform(project_to_ll, inter)
                center = cell.centroid
                center_ll = transform(project_to_ll, center)

                id_str = f"{prefix}_{idx:04d}"
                cells.append((inter_ll, Point(center_ll.x, center_ll.y), id_str))
                idx += 1

    return cells, idx


# ----------------------------
# Routes
# ----------------------------
@app.route("/", methods=["GET", "POST"])
def index():
    if request.method == "POST":
        file = request.files.get("kmlfile")
        grid_size_raw = request.form.get("grid_size")
        prefix = request.form.get("prefix") or "G"
        preview = request.form.get("preview") == "on"

        if not file or file.filename == "":
            flash("Please upload a KML/KMZ file.", "danger")
            return redirect(request.url)

        if not allowed_file(file.filename):
            flash("Invalid file type. Only .kml or .kmz allowed.", "danger")
            return redirect(request.url)

        try:
            grid_size = float(grid_size_raw)
        except:
            flash("Grid size must be a number.", "danger")
            return redirect(request.url)

        tmpdir = tempfile.mkdtemp(dir=app.config["UPLOAD_FOLDER"])
        filepath = os.path.join(tmpdir, file.filename)
        file.save(filepath)

        polys = parse_kml(filepath)

        if not polys:
            flash("No polygon found inside KML file.", "danger")
            return redirect(request.url)

        kml_out = simplekml.Kml()
        cur_idx = 1
        total_cells = 0

        for i, poly in enumerate(polys, start=1):
            cells, cur_idx = make_grid_for_polygon(poly, grid_size, cur_idx, prefix)
            total_cells += len(cells)

            folder = kml_out.newfolder(name=f"Polygon_{i}")

            for geom, center, cid in cells:
                pg = folder.newpolygon(
                    name=cid,
                    outerboundaryis=list(geom.exterior.coords)
                )
                pg.style.polystyle.fill = 0
                pg.style.linestyle.color = simplekml.Color.red

                pt = folder.newpoint(name=cid, coords=[(center.x, center.y)])
                pt.style.iconstyle.scale = 0

        out_path = os.path.join(tmpdir, "Grid.kml")
        kml_out.save(out_path)

        if preview:
            rel_path = os.path.basename(tmpdir) + "/Grid.kml"
            preview_url = url_for("download_file", path=rel_path)

            return render_template(
                "index.html",
                preview_url=preview_url,
                download_url=preview_url,
                total_cells=total_cells,
                grid_size=grid_size,
                prefix=prefix
            )

        return send_file(out_path, as_attachment=True)

    return render_template("index.html")


# ✅ Download endpoint
@app.route("/download")
def download_file():
    path = request.args.get("path")
    if not path:
        flash("Missing download path.", "danger")
        return redirect(url_for("index"))

    full = os.path.join(app.config["UPLOAD_FOLDER"], path)
    if not os.path.exists(full):
        flash("File not found.", "danger")
        return redirect(url_for("index"))

    return send_file(full, as_attachment=True)


# ✅ IMPORTANT: Do NOT add app.run() — Render uses gunicorn
# gunicorn runs this app as: gunicorn app:app


