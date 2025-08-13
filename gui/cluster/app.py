from flask import Flask, request, render_template, send_from_directory
from aptamerutils import SeqList
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
        header = request.form['header']
        header_tol = int(request.form['header_tol'])
        end = request.form['end']
        end_tol = int(request.form['end_tol'])
        eps = float(request.form['eps'])
        topk = int(request.form['topk'])
        displaykmerfeature = int(request.form['displaykmerfeature'])
        showLoop = request.form.get('showLoop') == 'on'
        featureminfrac = float(request.form['featureminfrac'])
        savetopk = int(request.form.get("savetopk"))
        fixedlength = int(request.form.get("fixedLengthInt"))
        tolerance = int(request.form.get("fixedLengthTol"))

        filepath = os.path.join(UPLOAD_FOLDER, file.filename)
        file.save(filepath)

        seq = SeqList().fromfastq(filepath)
        seq = seq.trimWithFuzzyPattern(header, end, fixedlength, tolerance, header_tol, end_tol)
        seq = seq.sortbyCount(topk = savetopk)
        clusters = seq.generateDistMap().createPosMap().getCluster(eps=eps)
        seq.getClustersLabeled(clusters)
        seq_sorted = seq.sortbyCount(topk=topk)
        seq_sorted = seq_sorted.sortbyCluster()

        outpath = os.path.join(RESULT_FOLDER, f'{os.path.splitext(file.filename)[0]}.pdf')
        seq_sorted.drawText(outpath,
                            displaykmerfeature=displaykmerfeature,
                            showLoop=showLoop,
                            featureminfrac=featureminfrac,
                            header=header,
                            end=end)

        return render_template('index.html', image_path = f"results/{os.path.splitext(file.filename)[0]}.pdf")

    return render_template('index.html', image_path=None)

@app.route('/results/<filename>')
def result_file(filename):
    return send_from_directory(RESULT_FOLDER, filename)

if __name__ == '__main__':
    app.run(debug=True)