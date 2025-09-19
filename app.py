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
def home():
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


@app.route("/blast", methods=["GET", "POST"])
def blast():
    blast_result = None
    submission = None

    if request.method == "POST":
        seq = request.form.get("sequence", "").strip()
        program = request.form.get("program", "blastn")
        db = request.form.get("db", "nt")

        if not seq:
            flash("Please provide a sequence to BLAST.", "warning")
            return redirect(url_for("blast"))

        # Show immediate submission summary
        submission = {"program": program, "db": db, "sequence": seq[:1000]}

        # ---------------------------
        # REMOTE BLAST: enabled here
        # ---------------------------
        # WARNING / PREPARATION:
        #   1) Set Entrez.email to your real email (NCBI policy).
        #   2) Use this sparingly and avoid heavy automation.
        #   3) Queries can take several seconds to return.
        #
        # Replace the placeholder email below with your real email address:
        Entrez.email = "pathangufran123786@gmail.com"   # <- REPLACE with your email before running!

        try:
            # qblast will send the query to NCBI and return an HTTP handle
            # format_type="XML" returns BLAST XML which we parse with NCBIXML
            handle = NCBIWWW.qblast(program, db, seq, format_type="XML", hitlist_size=10)
            blast_xml = handle.read()
            handle.close()

            # Parse BLAST XML: collect top hits/hsp info
            blast_records = NCBIXML.parse(io.StringIO(blast_xml))
            hits = []
            # NCBIXML.parse yields one or more records; we'll iterate and take alignments
            for record in blast_records:
                for alignment in record.alignments[:10]:  # top 10 alignments
                    for hsp in alignment.hsps[:1]:  # first HSP for each alignment
                        # compute percent identity (identities / align_length * 100) if available
                        ident = getattr(hsp, "identities", None)
                        alin_len = getattr(hsp, "align_length", None)
                        pct_id = None
                        if ident is not None and alin_len:
                            try:
                                pct_id = round((int(ident) / int(alin_len)) * 100, 2)
                            except Exception:
                                pct_id = None

                        hits.append({
                            "title": getattr(alignment, "title", ""),
                            "length": getattr(alignment, "length", None),
                            "score": getattr(hsp, "score", None),
                            "e_value": getattr(hsp, "expect", None),
                            "identity": ident,
                            "align_length": alin_len,
                            "pct_identity": pct_id
                        })
                # Since most qblast queries return a single record, break after first record
                break

            blast_result = {"hits": hits, "submission": submission}
            if not hits:
                blast_result["note"] = "No hits found."

        except Exception as e:
            # Catch network/NCBI errors and show friendly message
            flash(f"BLAST error: {e}", "danger")
            # Provide the submission summary so the user knows what was sent
            blast_result = {"error": str(e), "submission": submission}

    return render_template("blast.html", result=blast_result)


if __name__ == "__main__":
    app.run(debug=True, host="0.0.0.0", port=5000)
