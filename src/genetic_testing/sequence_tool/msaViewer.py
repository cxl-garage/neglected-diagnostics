import multiprocessing
import signal
import subprocess
import threading
from io import StringIO

import dash_bio as dashbio
from dash import Dash, html
from flask import Flask


def sequence_viewer(fasta_file, app):
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
    app.run(debug=False, use_reloader=False)


class MsaViewer(threading.Thread):
    def __init__(self, fasta_file):
        super().__init__()
        self.fasta_file = fasta_file
        self.app = Dash(__name__)
        self.viewerThread = threading.Thread(
            target=sequence_viewer,
            args=(
                self.fasta_file,
                self.app,
            ),
        )
        self.event = threading.Event()

    def updateData(self, fasta_file):
        self.fasta_file = fasta_file
        stringio = StringIO(fasta_file.getvalue().decode("utf-8"))
        data = stringio.read()
        self.app.layout = html.Div(
            [
                dashbio.AlignmentChart(
                    id="alignment-viewer",
                    data=data,
                ),
                html.Div(id="default-alignment-viewer-output"),
            ]
        )

    def stop(self):
        # Stop the Flask server
        self.viewerThread.join()

        self.event.set()

    def stopped(self):
        return self.event.is_set()
