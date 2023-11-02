FROM python:3.10
COPY ./src ./src
COPY ./requirements.txt ./src/requirements.txt
RUN pip install -r ./src/requirements.txt
EXPOSE 80
RUN mkdir /root/.streamlit
COPY deployment/credentials.toml /root/.streamlit/credentials.toml
COPY .streamlit/secrets.toml /root/.streamlit/secrets.toml
ENV PYTHONPATH ./src
ENTRYPOINT ["streamlit", "run"]
CMD ["src/app/Home.py", "--server.port", "80"]