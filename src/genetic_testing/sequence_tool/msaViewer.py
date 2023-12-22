import threading
from io import StringIO

import dash_bio as dashbio
from dash import Dash, html


def sequence_viewer(fasta_file):
    app = Dash(__name__)
    stringio = StringIO(fasta_file.getvalue().decode("utf-8"))
    data = stringio.read()
    print(type(data))
    app.layout = html.Div(
        [
            dashbio.AlignmentChart(
                id="alignment-viewer",
                data=data,
            ),
            html.Div(id="default-alignment-viewer-output"),
        ]
    )
    app.run_server(debug=False)


def start_sequence_viewer(fasta_file):
    t1 = threading.Thread(target=sequence_viewer, args=(fasta_file,))
    t1.start()
