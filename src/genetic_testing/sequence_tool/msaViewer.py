import threading

import dash_bio as dashbio
from dash import Dash, html


def sequence_viewer(sequences):
    app = Dash(__name__)
    data = sequences.decode("utf-8")
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


def start_sequence_viewer(sequences):
    t1 = threading.Thread(target=sequence_viewer, args=(sequences,))
    t1.start()
