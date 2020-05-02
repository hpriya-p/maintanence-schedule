FROM julia:latest
COPY src .
COPY instances .

# Add Julia libraries
ADD REQUIRE /.julia/v0.3/REQUIRE
RUN julia -e "Pkg.resolve()"