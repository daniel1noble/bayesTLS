#!/usr/bin/env python3
"""Strip the non-breaking space (U+00A0) Pandoc inserts between supplement
crossref prefixes ("Figure S", "Table S", "Equation S") and the number
that follows, inside a DOCX file (i.e. word/document.xml of the zip).

Used by the project Makefile after rendering supp.docx so the supplement
reads "Figure S1" rather than "Figure S<nbsp>1" — matching the standard
journal supplement format and the post-processed HTML and PDF outputs.

Usage:
    python3 scripts/strip-docx-nbsp.py path/to/supp.docx
"""

import os
import re
import shutil
import sys
import tempfile
import zipfile


def main(src: str) -> None:
    fd, tmp = tempfile.mkstemp(prefix="docx-strip-", suffix=".docx")
    os.close(fd)
    pattern = re.compile(rb"(Figure|Table|Equation) S\xc2\xa0")
    with zipfile.ZipFile(src, "r") as zin, \
            zipfile.ZipFile(tmp, "w", zipfile.ZIP_DEFLATED) as zout:
        for item in zin.infolist():
            data = zin.read(item.filename)
            if item.filename == "word/document.xml":
                data = pattern.sub(rb"\1 S", data)
            zout.writestr(item, data)
    shutil.move(tmp, src)


if __name__ == "__main__":
    if len(sys.argv) != 2:
        sys.exit("usage: strip-docx-nbsp.py <docx-file>")
    main(sys.argv[1])
