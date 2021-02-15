FROM python:3.6

RUN git clone https://github.com/churchmanlab/genewalk.git && \
    cd genewalk && \
    pip install . && \
    python -m genewalk.resources
