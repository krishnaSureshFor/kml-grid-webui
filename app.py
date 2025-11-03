# app.py — Full upgraded KML Grid WebUI
import os
import tempfile
import zipfile
import json
import io
import math
from datetime import datetime
from flask import Flask, render_template, request, send_file, redirect, url_for, flash, jsonify
from lxml import etree
from shapely.geometry import Polygon, Point, box, mapping, shape
from shapely.ops import transform, unary_union
import geopandas as gpd
import simplekml
from pyproj import Transformer, CRS
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
import shutil
import csv

# ----------------------------
# Config
# ----------------------------
UPLOAD_FOLDER = "uploads"
HISTORY_FILE = os.path.join(UPLOAD_FOLDER, "history.json")
os.makedirs(UPLOAD_FOLDER, exist_ok=True)

app = Flask(__name__)
app.secret_key = "secret123"
app.config["UPLOAD_FOLDER"] = UPLOAD_FOLDER
app.config["MAX_CONTENT_LENGTH"] = 50 * 1024 * 1024  # 50 MB

ALLOWED_EXT = {"kml", "kmz", "geojson", "json"}

# utility to persist job history
def load_history():
    if os.path.exists(HISTORY_FILE):
        try:
            return json.load(open(HISTORY_FILE))
        except Exception:
            return []
    return []

def save_history(history):
    json.dump(history, open(HISTORY_FILE, "w"), indent=2)

# ----------------------------
# Helpers
# ----------------------------
def allowed_file(filename):
    return "." in filename and filename.rsplit(".", 1)[1].lower() in ALLOWED_EXT

def mktempdir(prefix="job"):
    d = tempfile.mkdtemp(prefix=prefix + "_", dir=app.config["UPLOAD_FOLDER"])
    return d

def zip_shapefile(shp_basepath, out_zip):
    # shp_basepath is path without extension or path to .shp
    shp = shp_basepath if shp_basepath.endswith(".shp") else shp_basepath + ".shp"
    base = shp[:-4]
    exts = [".shp", ".shx", ".dbf", ".prj", ".cpg"]
    with zipfile.ZipFile(out_zip, "w", zipfile.ZIP_DEFLATED) as z:
        for e in exts:
            f = base + e
            if os.path.exists(f):
                z.write(f, os.path.basename(f))

def make_kmz(kml_path, kmz_path):
    with zipfile.ZipFile(kmz_path, "w", zipfile.ZIP_DEFLATED) as zf:
        zf.write(kml_path, arcname="doc.kml")

# ----------------------------
# KML parser (robust)
# ----------------------------
def parse_kml_or_geojson(filepath):
    """Returns list of shapely Polygons in lon/lat."""
    ext = filepath.lower().rsplit(".", 1)[-1]
    polygons = []
    if ext == "kmz":
        with zipfile.ZipFile(filepath, "r") as z:
            kml_files = [f for f in z.namelist() if f.lower().endswith(".kml")]
            if not kml_files:
                return polygons
            with z.open(kml_files[0]) as fh:
                xml = fh.read()
        root = etree.fromstring(xml)
        ns = {"kml": "http://www.opengis.net/kml/2.2"}
        coords_nodes = root.findall(".//kml:Polygon//kml:coordinates", ns)
        for node in coords_nodes:
            if node.text is None:
                continue
            raw = node.text.strip()
            parts = raw.replace("\n", " ").replace("\t", " ").split()
            ring = []
            for p in parts:
                vals = p.split(",")
                if len(vals) >= 2:
                    try:
                        lon = float(vals[0]); lat = float(vals[1])
                        ring.append((lon, lat))
                    except: pass
            if len(ring) >= 3:
                polygons.append(Polygon(ring))
        return polygons
    elif ext in ("kml",):
        with open(filepath, "rb") as fh:
            xml = fh.read()
        root = etree.fromstring(xml)
        ns = {"kml": "http://www.opengis.net/kml/2.2"}
        coords_nodes = root.findall(".//kml:Polygon//kml:coordinates", ns)
        for node in coords_nodes:
            if node.text is None: continue
            raw = node.text.strip()
            parts = raw.replace("\n", " ").replace("\t", " ").split()
            ring = []
            for p in parts:
                vals = p.split(",")
                if len(vals) >= 2:
                    try:
                        lon = float(vals[0]); lat = float(vals[1])
                        ring.append((lon, lat))
                    except: pass
            if len(ring) >= 3:
                polygons.append(Polygon(ring))
        return polygons
    else:
        # geojson or json
        try:
            gdf = gpd.read_file(filepath)
            for geom in gdf.geometry:
                if geom is None: continue
                if geom.geom_type == "Polygon":
                    polygons.append(geom)
                elif geom.geom_type == "MultiPolygon":
                    for g in geom.geoms:
                        polygons.append(g)
            return polygons
        except Exception:
            return polygons

# ----------------------------
# UTM helpers
# ----------------------------
def utm_crs_for_lonlat(lon, lat):
    # returns EPSG code for UTM zone (north/south)
    zone = int((lon + 180) / 6) + 1
    if lat >= 0:
        epsg = 32600 + zone
    else:
        epsg = 32700 + zone
    return CRS.from_epsg(epsg)

# ----------------------------
# Grid generation with UTM per polygon
# ----------------------------
def make_grid_for_polygon_utm(poly_ll, grid_m, id_start, prefix="G"):
    """
    For accurate square grids, detect suitable UTM zone based on polygon centroid.
    Returns list of (geom_lonlat, area_m2, center_lonlat, id_str)
    """
    centroid = poly_ll.centroid
    lon, lat = centroid.x, centroid.y
    utm_crs = utm_crs_for_lonlat(lon, lat)
    crs_from = CRS.from_epsg(4326)
    project_to_m = Transformer.from_crs(crs_from, utm_crs, always_xy=True).transform
    project_to_ll = Transformer.from_crs(utm_crs, crs_from, always_xy=True).transform

    poly_m = transform(project_to_m, poly_ll)
    minx, miny, maxx, maxy = poly_m.bounds

    cells = []
    idx = id_start

    nx = max(1, int((maxx - minx) / grid_m) + 1)
    ny = max(1, int((maxy - miny) / grid_m) + 1)

    for i in range(nx):
        for j in range(ny):
            x0 = minx + i * grid_m
            y0 = miny + j * grid_m
            cell = box(x0, y0, x0 + grid_m, y0 + grid_m)
            if poly_m.intersects(cell):
                inter_m = poly_m.intersection(cell)
                if inter_m.is_empty:
                    continue
                # area in m2
                try:
                    area_m2 = float(inter_m.area)
                except:
                    area_m2 = 0.0
                inter_ll = transform(project_to_ll, inter_m)
                center_m = cell.centroid
                center_ll = transform(project_to_ll, center_m)
                id_str = f"{prefix}_{idx:05d}"
                cells.append((inter_ll, area_m2, Point(center_ll.x, center_ll.y), id_str))
                idx += 1

    return cells, idx

# ----------------------------
# Export helpers: KML, KMZ, GeoJSON, Shapefile, GPKG, CSV, PDF
# ----------------------------
def save_kml_from_cells(cells, out_kml_path, description="Grid"):
    kml_out = simplekml.Kml()
    folder = kml_out.newfolder(name="Grids")
    for geom, area_m2, center, gid in cells:
        coords = list(geom.exterior.coords)
        pol = folder.newpolygon(name=gid, outerboundaryis=coords)
        pol.style.polystyle.fill = 0
        pol.style.linestyle.width = 1
    kml_out.document.description = description
    kml_out.save(out_kml_path)

def save_geojson_from_cells(cells, out_geojson_path):
    features = []
    geoms = []
    ids = []
    areas = []
    for geom, area_m2, center, gid in cells:
        geoms.append(mapping(geom))
        ids.append(gid)
        areas.append(area_m2/10000.0)
    gdf = gpd.GeoDataFrame({"Grid_ID": ids, "Area_ha": areas, "geometry": [shape(g) for g in geoms]}, crs="EPSG:4326")
    gdf.to_file(out_geojson_path, driver="GeoJSON")

def save_shapefile_from_cells(cells, out_shp_base):
    geoms = []
    ids = []
    areas = []
    for geom, area_m2, center, gid in cells:
        geoms.append(geom)
        ids.append(gid)
        areas.append(area_m2/10000.0)
    gdf = gpd.GeoDataFrame({"Grid_ID": ids, "Area_ha": areas, "geometry": geoms}, crs="EPSG:4326")
    # write shapefile at out_shp_base + .shp
    gdf.to_file(out_shp_base + ".shp", driver="ESRI Shapefile")

def save_gpkg_from_cells(cells, out_gpkg):
    geoms = []
    ids = []
    areas = []
    for geom, area_m2, center, gid in cells:
        geoms.append(geom)
        ids.append(gid)
        areas.append(area_m2/10000.0)
    gdf = gpd.GeoDataFrame({"Grid_ID": ids, "Area_ha": areas, "geometry": geoms}, crs="EPSG:4326")
    gdf.to_file(out_gpkg, driver="GPKG")

def save_csv_from_cells(cells, out_csv):
    with open(out_csv, "w", newline='', encoding='utf-8') as fh:
        writer = csv.writer(fh)
        writer.writerow(["Grid_ID", "Area_ha"])
        for geom, area_m2, center, gid in cells:
            writer.writerow([gid, round(area_m2/10000.0, 4)])

def save_pdf_map(cells, out_pdf, title="Grid Map"):
    # Simple plot: polygons with labels (optionally)
    fig, ax = plt.subplots(figsize=(10,10))
    gdf = gpd.GeoDataFrame([{"geometry": geom, "Grid_ID": gid} for geom, area_m2, center, gid in cells], crs="EPSG:4326")
    # project to mercator for basemap compatibility
    try:
        gdf_m = gdf.to_crs(epsg=3857)
        gdf_m.boundary.plot(ax=ax, linewidth=0.5)
    except Exception:
        gdf.boundary.plot(ax=ax, linewidth=0.5)
    ax.set_title(title)
    ax.set_axis_off()
    pp = PdfPages(out_pdf)
    pp.savefig(fig, bbox_inches='tight')
    pp.close()
    plt.close(fig)

# ----------------------------
# Main route
# ----------------------------
@app.route("/", methods=["GET", "POST"])
def index():
    if request.method == "POST":
        # form values
        use_drawn = request.form.get("use_drawn") == "on"
        merge_polys = request.form.get("merge_polys") == "on"
        show_labels = request.form.get("show_labels") == "on"
        export_geojson = request.form.get("export_geojson") == "on"
        export_shp = request.form.get("export_shp") == "on"
        export_gpkg = request.form.get("export_gpkg") == "on"
        export_pdf = request.form.get("export_pdf") == "on"

        # uploaded file or drawn geojson (drawn geojson comes in request.files as 'drawn')
        file = request.files.get("kmlfile")
        drawn_geojson = request.form.get("drawn_geojson")

        # grid size, prefix
        grid_size_raw = request.form.get("grid_size")
        prefix = request.form.get("prefix") or "G"
        preview = request.form.get("preview") == "on"

        # validations
        if not file and not drawn_geojson:
            flash("Upload KML/KMZ/GeoJSON or draw polygon on map.", "danger")
            return redirect(request.url)

        try:
            grid_size = float(grid_size_raw)
            if grid_size <= 0: raise ValueError()
        except:
            flash("Grid size must be positive number (meters).", "danger")
            return redirect(request.url)

        # prepare tmpdir
        tmpdir = mktempdir("job")
        cells_all = []  # flattened cells for output
        cur_idx = 1

        # 1) Read polygons
        polys = []
        if drawn_geojson:
            try:
                gj = json.loads(drawn_geojson)
                # features maybe present
                if "features" in gj:
                    for feat in gj["features"]:
                        geom = shape(feat["geometry"])
                        if geom.geom_type == "Polygon":
                            polys.append(geom)
                        elif geom.geom_type == "MultiPolygon":
                            for g in geom.geoms: polys.append(g)
                else:
                    geom = shape(gj)
                    if geom.geom_type == "Polygon":
                        polys.append(geom)
                    elif geom.geom_type == "MultiPolygon":
                        for g in geom.geoms: polys.append(g)
            except Exception as e:
                flash("Failed to parse drawn GeoJSON: " + str(e), "danger")
                return redirect(request.url)

        if file:
            fname = file.filename
            fpath = os.path.join(tmpdir, fname)
            file.save(fpath)
            parsed = parse_kml_or_geojson(fpath)
            if parsed:
                polys.extend(parsed)

        if not polys:
            flash("No polygons available to grid.", "danger")
            return redirect(request.url)

        # optionally merge polygons
        if merge_polys:
            try:
                unioned = unary_union(polys)
                merged = []
                if unioned.geom_type == "Polygon":
                    merged = [unioned]
                elif unioned.geom_type == "MultiPolygon":
                    merged = list(unioned.geoms)
                else:
                    merged = polys
                polys = merged
            except Exception:
                pass

        # generate grids per polygon using UTM detection
        all_cells = []
        for poly in polys:
            cells, cur_idx = make_grid_for_polygon_utm(poly, grid_size, cur_idx, prefix)
            all_cells.extend(cells)

        if not all_cells:
            flash("Grid generation produced zero cells (maybe grid too large).", "danger")
            return redirect(request.url)

        # Save outputs
        # KML
        out_kml = os.path.join(tmpdir, "Grid.kml")
        save_kml_from_cells(all_cells, out_kml, description=f"UI Developed by Krishnasureshkumar FG")

        # KMZ
        out_kmz = os.path.join(tmpdir, "Grid.kmz")
        try:
            make_kmz(out_kml, out_kmz)
        except:
            out_kmz = None

        # GeoJSON
        out_geojson = os.path.join(tmpdir, "Grid.geojson")
        try:
            save_geojson_from_cells(all_cells, out_geojson)
        except:
            out_geojson = None

        # Shapefile (zipped)
        out_shp_base = os.path.join(tmpdir, "Grid")
        try:
            save_shapefile_from_cells(all_cells, out_shp_base)
            out_shp_zip = os.path.join(tmpdir, "Grid_shapefile.zip")
            zip_shapefile(out_shp_base, out_shp_zip)
        except:
            out_shp_zip = None

        # GPKG
        out_gpkg = os.path.join(tmpdir, "Grid.gpkg")
        try:
            save_gpkg_from_cells(all_cells, out_gpkg)
        except:
            out_gpkg = None

        # CSV attribute
        out_csv = os.path.join(tmpdir, "Grid_Table.csv")
        try:
            save_csv_from_cells(all_cells, out_csv)
        except:
            out_csv = None

        # PDF map
        out_pdf = os.path.join(tmpdir, "Grid_Map.pdf")
        try:
            if export_pdf:
                save_pdf_map(all_cells, out_pdf, title="Grid Export")
        except:
            out_pdf = None

        # record job history
        history = load_history()
        job = {
            "id": os.path.basename(tmpdir),
            "time": datetime.utcnow().isoformat(),
            "cells": len(all_cells),
            "kml": "Grid.kml" if os.path.exists(out_kml) else None,
            "kmz": "Grid.kmz" if os.path.exists(out_kmz) else None,
            "geojson": "Grid.geojson" if os.path.exists(out_geojson) else None,
            "shp": "Grid_shapefile.zip" if os.path.exists(os.path.join(tmpdir, "Grid_shapefile.zip")) else None,
            "gpkg": "Grid.gpkg" if os.path.exists(out_gpkg) else None,
            "csv": "Grid_Table.csv" if os.path.exists(out_csv) else None,
            "pdf": "Grid_Map.pdf" if os.path.exists(out_pdf) else None
        }
        history.insert(0, job)
        # keep only last 20
        history = history[:20]
        save_history(history)

        # build download urls (relative)
        def r(p):
            return url_for("download_file", path=os.path.basename(tmpdir) + "/" + p) if p else None

        preview_url = r("Grid.kml")
        return render_template(
            "index.html",
            preview_url=preview_url,
            download_url=r("Grid.kml"),
            download_kmz=r("Grid.kmz"),
            download_geojson=r("Grid.geojson"),
            download_shp=r("Grid_shapefile.zip"),
            download_gpkg=r("Grid.gpkg"),
            download_csv=r("Grid_Table.csv"),
            download_pdf=r("Grid_Map.pdf"),
            total_cells=len(all_cells)
        )

    # GET
    history = load_history()
    return render_template("index.html", history=history)


# ----------------------------
# Download endpoint
# ----------------------------
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

# Note: Do not use app.run() on Render; use gunicorn app:app
