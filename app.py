import os
import io
import threading
import uuid

from flask import Flask, render_template, request, redirect, url_for, flash, jsonify
from Bio.Seq import Seq
from Bio.Data import CodonTable
from Bio.Blast import NCBIWWW, NCBIXML
from Bio import Entrez

# Read email & secret key from environment (required by NCBI + for security)
ENTREZ_EMAIL = os.environ.get("ENTREZ_EMAIL", "pathangufran123786@gmail.com")
SECRET_KEY = os.environ.get("SECRET_KEY", "replace-this-with-a-secure-random-key")

# Set Entrez email globally (NCBI requires a valid email)
Entrez.email = ENTREZ_EMAIL

app = Flask(__name__)
app.secret_key = SECRET_KEY  # Flask sessions / flash messages


 # Simple AA -> codon mapping using the first codon for each amino acid
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


import threading
import uuid
from flask import jsonify

# Simple in-memory job store
JOBS = {}  # job_id -> {"status": "running"|"done"|"error", "result": "<xml/text>"}

def run_blast(job_id, seq, program="blastn", database="nt"):
    """
    Run BLAST in a background thread and save result in JOBS.
    """
    try:
        # Use the Entrez email you already configured earlier
        from Bio import Entrez
        Entrez.email = os.environ.get("ENTREZ_EMAIL", "pathangufran123786@gmail.com")

        handle = NCBIWWW.qblast(program, database, seq)
        result_text = handle.read()
        handle.close()

        JOBS[job_id]["status"] = "done"
        JOBS[job_id]["result"] = result_text
    except Exception as e:
        JOBS[job_id]["status"] = "error"
        JOBS[job_id]["result"] = str(e)

@app.route("/blast", methods=["GET", "POST"])
def blast():
    if request.method == "GET":
        return render_template("blast.html")

    seq = request.form.get("sequence", "").strip()
    if not seq:
        return jsonify({"error": "no sequence provided"}), 400

    job_id = str(uuid.uuid4())
    JOBS[job_id] = {"status": "running", "result": None}

    # Start BLAST in a background thread
    t = threading.Thread(target=run_blast, args=(job_id, seq))
    t.daemon = True
    t.start()

    # Return job id immediately (JSON)
    return jsonify({"job_id": job_id}), 202

@app.route("/blast_status/<job_id>", methods=["GET"])
def blast_status(job_id):
    job = JOBS.get(job_id)
    if not job:
        return jsonify({"error": "job not found"}), 404
    return jsonify(job)


if __name__ == "__main__":
    app.run(debug=True, host="0.0.0.0", port=5000)
