FROM python:3.10
COPY dist/*.whl .
RUN pip install *.whl
EXPOSE 8501
COPY ./src ./src
RUN mkdir /root/.streamlit
COPY deployment/credentials.toml /root/.streamlit/credentials.toml
COPY .streamlit/secrets.toml /root/.streamlit/secrets.toml
ENTRYPOINT ["streamlit", "run"]
CMD ["src/app/main.py"]