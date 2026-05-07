# CIGALE galaxy/AGN SED-fitting mini-project package

Contents:

- `sample_targets.csv` — four-galaxy sample for the assignment.
- `make_cigale_photometry_from_dustpedia.py` — script that queries the DustPedia VizieR table and writes a CIGALE input file.
- `pcigale_AGN.ini` — starter CIGALE configuration using an AGN torus component.
- `pcigale_noAGN.ini` — same idea, but without AGN; useful for comparison.
- `assignment_writeup.txt` — student-facing assignment text.
- `estimate_optical_agn_fraction.py` — optional helper script to estimate AGN/total model fraction at optical wavelengths from best-model SED files.

Recommended workflow:

1. Install CIGALE and helper packages in a fresh conda environment:

   conda create -n cigale-class python=3.12 -y
   conda activate cigale-class
   pip install pcigale astropy pyvo pandas numpy matplotlib

2. Generate the photometry table:

   python make_cigale_photometry_from_dustpedia.py

   This should produce `cigale_dustpedia_input.csv`.

3. Test whether CIGALE sees the filters:

   pcigale-filters list | grep -i galex
   pcigale-filters list | grep -i wise

4. Run the no-AGN and AGN fits in separate directories:

   mkdir run_noAGN run_AGN
   cp cigale_dustpedia_input.csv pcigale_noAGN.ini run_noAGN/
   cp cigale_dustpedia_input.csv pcigale_AGN.ini run_AGN/
   cd run_noAGN
   cp pcigale_noAGN.ini pcigale.ini
   pcigale check
   pcigale run
   pcigale-plots sed
   cd ../run_AGN
   cp pcigale_AGN.ini pcigale.ini
   pcigale check
   pcigale run
   pcigale-plots sed

5. Compare outputs in `out/results.fits` and the SED plots in `out/`.

Notes:

- DustPedia photometry is aperture-matched whole-galaxy photometry. That is good pedagogically, but it means the AGN is diluted by host galaxy light for some objects.
- NGC 4151 is variable, so non-simultaneous broadband photometry should not be over-interpreted.
- If the installed CIGALE version uses revised filter names, edit the column names in the generated CSV and the `bands =` line in the `.ini` files. CIGALE 2025 introduced a more systematic filter naming scheme.
