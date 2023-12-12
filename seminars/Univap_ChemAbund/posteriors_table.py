from pathlib import Path
from delete.data_printing import PdfPrinter, latex_labels
from pylatex import NoEscape

atomic_data_dict = {'T_low': r'$Normal(\mu=15,000K,\sigma=5000K)$',
                    'T_high': r'$Normal(\mu=15,000K,\sigma=5000K)$',
                    'n_e': r'$HalfCauchy(\mu=2.0,\sigma=0)$',
                    'cHbeta': r'$HalfCauchy(\mu=2.0,\sigma=0)$',
                    'log(X_i+)': r'$Normal(\mu=5,\sigma=5)$',
                    'log(He1r)': r'$Normal(\mu=0,\sigma=3)$',
                    'log(He2r)': r'$Normal(\mu=0,\sigma=3)$',
                    'Teff': r'$Uniform(min=30,000K, max=90,000K)$',
                    'logU': r'$Uniform(min=-4.0, max=-1.0K)$'}

tables_folder = Path('/home/vital/Dropbox/Astrophysics/Seminars/UniVapo 2021/')

pdf = PdfPrinter()

pdf.create_pdfDoc(pdf_type='table')
pdf.pdfDoc.append(NoEscape('\definecolor{background}{rgb}{0.169, 0.169, 0.169}'))
pdf.pdfDoc.append(NoEscape('\definecolor{foreground}{rgb}{0.702, 0.780, 0.847}'))
pdf.pdfDoc.append(NoEscape(r'\arrayrulecolor{foreground}'))

pdf.pdf_insert_table(['Parameter', 'Prior distribution'],
                     addfinalLine=True, color_font='foreground', color_background='background')

for ion, params in atomic_data_dict.items():
    ion_label = latex_labels[ion]
    row = [ion_label, params]
    pdf.addTableRow(row, last_row=False, color_font='foreground', color_background='background')

pdf.table.add_hline()
pdf.generate_pdf(tables_folder/'posteriors_table', clean_tex=True)

