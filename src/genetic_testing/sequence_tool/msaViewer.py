import configparser
import os
import shutil
import subprocess
import threading
from io import StringIO

import dash_bio as dashbio
from dash import Dash, html


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


def make_archive(source, destination):
    base = os.path.basename(destination)
    name = base.split(".")[0]
    format = base.split(".")[1]
    archive_from = os.path.dirname(source)
    archive_to = os.path.basename(source.strip(os.sep))
    shutil.make_archive(name, format, archive_from, archive_to)
    shutil.move("%s.%s" % (name, format), destination)


def msa_cleaner(file, config):
    filename = file.name.split(".")[0]
    zipFileName = filename + "_cleaned.zip"
    tempFile = "temp.fasta"
    tempConfig = "temp.ini"
    outDirName = "msa_cleaner_out"
    outdir = os.path.join(os.getcwd(), outDirName)
    if not os.path.exists(outdir):
        os.mkdir(outdir)
    outFile = os.path.join(outdir, filename)

    with open(tempFile, "wb") as f:
        f.write(file.read())

    configPath = os.path.join(os.getcwd(), "static/msa_cleaner_default.ini")
    shutil.copyfile(configPath, tempConfig)
    actConfig = configparser.ConfigParser()
    actConfig.read(tempConfig)
    for section in config.sections():
        for option in config.options(section):
            if option in actConfig[section]:
                actConfig[section][option] = config[section][option]

    with open(tempConfig, "w") as f:
        actConfig.write(f)

    proc = subprocess.Popen(
        [
            "CIAlign",
            "--inifile",
            tempConfig,
            "--infile",
            tempFile,
            "--outfile_stem",
            outFile,
        ],
        stdout=subprocess.PIPE,
        stderr=subprocess.DEVNULL,
    )
    print(proc.stdout.read())
    os.remove(tempFile)
    os.remove(tempConfig)
    zipPath = os.path.join(os.getcwd(), zipFileName)
    make_archive(outdir, zipPath)
    shutil.rmtree(outdir)
    return zipPath


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

    def run(self):
        self.viewerThread.start()

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
