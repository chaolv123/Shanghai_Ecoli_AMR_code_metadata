# Public Release Checklist

Before uploading or making the GitHub repository public:

- Confirm that farm and district metadata are acceptable for public release.
- Confirm that figures, supplementary tables, manuscript files and generated
  output workbooks are not present.
- Confirm that no raw sequencing reads are present.
- Confirm that `LICENSE` and `DATA_LICENSE.md` are present.
- Add final raw-read accession or repository links to `README.md` when
  available.
- Run:

```bash
pip install -r requirements.txt
python code/run_all.py
python -m unittest discover -s tests -v
```
