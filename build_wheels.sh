# need to create a requirements file before starting docker
# pip install pipreqs
# pipreqs --force . # --ignore=(ls | grep egg.info) . # creating requirements .txt
# sed -i '' '/\.egg==info/d' requirements.txt

# docker pull quay.io/pypa/manylinux1_x86_64
# # or docker pull quay.io/pypa/manylinux1_i686
# docker run -it -v (pwd):/io quay.io/pypa/manylinux1_x86_64

# yum install zlib-devel

for PYBIN in /opt/python/*[5-7]*/bin; do
    "${PYBIN}/pip" install cython numpy pysam # install these requirements first
    "${PYBIN}/pip" install -r /io/requirements.txt
    "${PYBIN}/pip" wheel /io/ -w wheelhouse/
done

for whl in wheelhouse/*.whl; do
    auditwheel repair "$whl" -w /io/wheelhouse/
done
