# Should be run from the root of the project
python -m build .
docker build -t ndiag --platform=linux/arm64 -f deployment/Dockerfile .

# docker run --rm -p 127.0.0.1:8501:8501 ndiag
