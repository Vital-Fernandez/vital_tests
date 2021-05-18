from pathlib import Path
from src.specsiser.data_printing import PdfPrinter, latex_labels
from pylatex import NoEscape

atomic_data_dict = {'O2': ['Pradhan et al. (2006); Tayal (2007)', 'Zeippen (1982); Wiese et al. (1996)'],
                    'O3': ['Aggarwal \& Keenan (2000)',            'Storey \& Zeippen (2000); Wiese et al. (1996)'],
                    'N2': ['Tayal (2011)',                        'Wiese et al. (1996); Galav\'is et al. (1997)'],
                    'S2': ['Tayal \& Zatsarinny (2010)',           'Podobedova et al. (2009)'],
                    'S3': ['Hudson et al. (2012)',                'Podobedova et al. (2009)'],
                    'S4': ['Tayal (2000)',                        'Dufton et al. (1982); Johnson et al. (1986)'],
                    'Ar3': ['Galav\'is et al. (1995)',              'Kaufman \& Sugar (1986); Galav\'is et al. (1995)'],
                    'Ar4': ['Ramsbottom \& Bell (1997)',           'Mendoza \& Zeippen (1982)'],
                    'Cl3': ['Butler \& Zeippen (1989)',            'Mendoza (1983)'],
                    'Fe3': ['Zhang (1996)',                       'Quinet (1996), Johansson et al. (2000)'],
                    'Ni3': ['Bautista (2001)',                    'Bautista (2001)']}

tables_folder = Path('/home/vital/Dropbox/Astrophysics/Seminars/UniVapo 2021/')

table_headers = ['SDSS ID', 'Label', 'R.A', 'DEC', 'z', '$Grisms^{1}$', '$T_{exp}$', 'Seeing', r'Standard stars$^{2}$']

pdf = PdfPrinter()

pdf.create_pdfDoc(pdf_type='table')
pdf.pdfDoc.append(NoEscape('\definecolor{background}{rgb}{0.169, 0.169, 0.169}'))
pdf.pdfDoc.append(NoEscape('\definecolor{foreground}{rgb}{0.702, 0.780, 0.847}'))
pdf.pdfDoc.append(NoEscape(r'\arrayrulecolor{foreground}'))

pdf.pdf_insert_table(['Ion', 'Collision Strengths', 'Transition probabilities'],
                     addfinalLine=True, color_font='foreground', color_background='background')

for ion, params in atomic_data_dict.items():
    ion_label = latex_labels[ion]
    row = [ion_label] + params
    pdf.addTableRow(row, last_row=False, color_font='foreground', color_background='background')

pdf.table.add_hline()
pdf.generate_pdf(tables_folder/'atomic_data', clean_tex=True)