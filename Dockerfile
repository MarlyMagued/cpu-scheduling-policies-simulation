FROM gcc:12
WORKDIR /usr/src/app
COPY . .
RUN make
CMD ["./lab6"]
