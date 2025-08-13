from flask import Flask, request, render_template, send_from_directory
from aptamerutils import SeqList, Find
import os

app = Flask(__name__)
UPLOAD_FOLDER = os.path.join(app.root_path, 'uploads')
RESULT_FOLDER = os.path.join(app.root_path, 'static', 'results')
app.config['UPLOAD_FOLDER'] = UPLOAD_FOLDER

os.makedirs(UPLOAD_FOLDER, exist_ok=True)
os.makedirs(RESULT_FOLDER, exist_ok=True)

@app.route('/', methods=['GET', 'POST'])
def index():
    if request.method == 'POST':
        file = request.files['fastq']
        trim = request.form.get("trim") == "on"
        if trim:
            header = request.form['header']
            end = request.form['end']
            headertol = int(request.form["header_tol"])
            endtol = int(request.form["end_tol"])
            expectedLength = int(request.form["fixedLengthInt"])
            lengthtol = int(request.form["fixedLengthTol"])
        pattern = request.form["customfeature"]

        filepath = os.path.join(UPLOAD_FOLDER, file.filename)
        file.save(filepath)
        seq = SeqList().fromfastq(filepath)
        if trim:
            seq = seq.trimWithFuzzyPattern(header, end, expectedLength, lengthtol, headertol, endtol)
        matchseq = seq.findCustomSeqCombination(eval(pattern))
        result_file = os.path.join(RESULT_FOLDER, f'{os.path.splitext(file.filename)[0]}.txt')
        matchseq.saveLines(result_file)
        
        return render_template('index.html', result_file = f'results/{os.path.splitext(file.filename)[0]}.txt')

    return render_template('index.html', result_file=None)

@app.route('/results/<filename>')
def result_file(filename):
    return send_from_directory(RESULT_FOLDER, filename)

if __name__ == '__main__':
    app.run(debug=True)