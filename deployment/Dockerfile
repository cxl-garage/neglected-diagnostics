FROM python:3.10
# Install dependencies through the wheel. Installing from source (pip install .) 
# fails because there's no git in the container so fetching the version from SCM fails.
COPY dist/*.whl .
RUN pip install *.whl
EXPOSE 80
COPY ./src ./src
RUN mkdir /root/.streamlit
COPY deployment/credentials.toml /root/.streamlit/credentials.toml
COPY .streamlit/secrets.toml /root/.streamlit/secrets.toml
ENTRYPOINT ["streamlit", "run"]
CMD ["src/app/Home.py", "--server.port", "80"]