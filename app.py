from flask import Flask, render_template, request, redirect, url_for, flash
from Bio.Seq import Seq
from Bio.Data import CodonTable
# The following imports are used for remote BLAST
from Bio.Blast import NCBIWWW, NCBIXML
from Bio import Entrez
import io

app = Flask(__name__)
app.secret_key = "replace-this-with-a-secure-random-key"  # change in production

#
# Simple AA -> codon mapping using the first codon for each amino acid
#
standard_table = CodonTable.unambiguous_dna_by_name["Standard"]
aa_to_codons = {}
for codon, aa in standard_table.forward_table.items():
    aa_to_codons.setdefault(aa, []).append(codon)
stop_codons = standard_table.stop_codons  # list of stop codons


@app.route("/")
def index():
    return render_template("index.html")


@app.route("/nt_to_aa", methods=["GET", "POST"])
def nt_to_aa():
    result = None
    input_seq = ""
    frame = 1
    translate_table = 1
    keep_stops = False

    if request.method == "POST":
        input_seq = request.form.get("sequence", "").strip()
        frame = int(request.form.get("frame", "1"))
        translate_table = int(request.form.get("table", "1"))
        keep_stops = bool(request.form.get("keep_stops", None))

        seq = "".join(input_seq.split()).upper()
        seq = seq.replace("U", "T")  # normalize RNA -> DNA

        if not seq:
            flash("Please paste a nucleotide sequence.", "warning")
            return redirect(url_for("nt_to_aa"))

        allowed = set("ATGCN")
        if not set(seq).issubset(allowed):
            flash("Sequence contains invalid characters. Only A/T/G/C/U/N allowed.", "danger")
            return render_template("nt_to_aa.html", result=None, input_seq=input_seq, frame=frame, table=translate_table, keep_stops=keep_stops)

        try:
            dna = Seq(seq)
            start = max(0, frame - 1)
            protein = dna[start:].translate(table=translate_table, to_stop=(not keep_stops))
            result = str(protein)
        except Exception as e:
            flash(f"Error translating sequence: {e}", "danger")

    return render_template("nt_to_aa.html", result=result, input_seq=input_seq, frame=frame, table=translate_table, keep_stops=keep_stops)


@app.route("/aa_to_nt", methods=["GET", "POST"])
def aa_to_nt():
    result = None
    input_seq = ""
    if request.method == "POST":
        input_seq = request.form.get("sequence", "").strip().upper()
        if not input_seq:
            flash("Please provide an amino acid sequence.", "warning")
            return redirect(url_for("aa_to_nt"))

        try:
            nts = []
            for aa in input_seq:
                if aa == "*":
                    codon = stop_codons[0] if stop_codons else "TAA"
                elif aa in aa_to_codons:
                    codon = aa_to_codons[aa][0]
                else:
                    codon = "NNN"
                nts.append(codon)
            result = "".join(nts)
        except Exception as e:
            flash(f"Error reverse-translating sequence: {e}", "danger")

    return render_template("aa_to_nt.html", result=result, input_seq=input_seq)


@app.route("/gc_content", methods=["GET", "POST"])
def gc_content():
    result = None
    input_seq = ""
    if request.method == "POST":
        input_seq = request.form.get("sequence", "").strip()
        seq = "".join(input_seq.split()).upper()
        seq = seq.replace("U", "T")
        if not seq:
            flash("Please provide a nucleotide sequence.", "warning")
            return redirect(url_for("gc_content"))

        valid = [c for c in seq if c in "ATGC"]
        total = len(valid)
        if total == 0:
            flash("No valid nucleotide characters found (A/T/G/C).", "warning")
            return redirect(url_for("gc_content"))

        gc_count = sum(1 for c in valid if c in "GC")
        gc_percent = (gc_count / total) * 100
        at_percent = 100 - gc_percent
        result = {
            "length": total,
            "gc_count": gc_count,
            "gc_percent": round(gc_percent, 3),
            "at_percent": round(at_percent, 3)
        }

    return render_template("gc_content.html", result=result, input_seq=input_seq)


@app.route('/blast', methods=['GET', 'POST'])
def blast():
    if request.method == 'GET':
        return render_template('blast.html')

    # POST request: user submitted a sequence
    seq = request.form.get('sequence', '').strip()
    if not seq:
        return render_template('blast.html', error="Please provide a sequence.")

    try:
        # Try running BLAST against a smaller DB to reduce timeouts
        result_handle = NCBIWWW.qblast("blastn", "refseq_rna", seq, format_type="XML", hitlist_size=5)
        xml_result = result_handle.read()
        result_handle.close()

        # Send raw result to a simple results page
        return render_template('blast_results.html', raw_xml=xml_result)

    except Exception as e:
        # Print full traceback to logs
        traceback.print_exc(file=sys.stdout)

        # Show user-friendly message
        return render_template(
            'blast.html',
            error="BLAST request failed. Possible reasons: NCBI timeout, network issue, or invalid input."
        )


if __name__ == "__main__":
    app.run(debug=True, host="0.0.0.0", port=5000)
