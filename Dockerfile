ARG PGX_DOCKER_TAG=latest
FROM 820411415250.dkr.ecr.us-east-1.amazonaws.com/helix-py-app-pgx:${PGX_DOCKER_TAG}

WORKDIR /pgx-verification
ENV PYTHONPATH=/src:/pgx-verification
ENV TZ=America/Los_Angeles

RUN apt-get update && apt-get -y upgrade

COPY . .

RUN pip3 install --upgrade pip && \
    pip3 install --no-cache-dir --target /src -r requirements.txt \
    pip3 install awscli --upgrade

# Override the pgx lambda entrypoint to allow for interactive shell
ENTRYPOINT ["/usr/bin/env"]
CMD ["/bin/bash"]
