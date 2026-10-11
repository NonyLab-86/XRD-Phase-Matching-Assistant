# XRD Phase-Matching Assistant — HTML Edition

A static browser app. Upload up to 5 experimental XRD files and 20 reference files; processing occurs locally in the browser.

## Publish with GitHub Pages
1. Create a public GitHub repository.
2. Upload `index.html` to the repository root.
3. Go to **Settings → Pages**.
4. Choose **Deploy from a branch**, branch `main`, folder `/(root)`, then Save.
5. Open the URL shown in Pages settings and test uploads.

## Supported files
Experimental: XY, TXT, DAT, CSV, XLSX, XLS. References: RTF, TXT, DAT, XY, CSV, XLSX, XLS. XLSX support uses SheetJS from a CDN, so an internet connection is required.

## Limitations
This version does not perform ALS baseline correction or robust PDF extraction/OCR. Export reference PDF peak tables to TXT/RTF first. Peak matching is a screening tool, not definitive phase identification. Check reference quality, diagnostic peaks, overlap, composition/EDS, and consider profile fitting/Rietveld refinement. Do not redistribute licensed reference database files without permission.
