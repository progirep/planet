set -e
unset -v latest
for file in "snapshots"/_iter_*.caffemodel; do
  [[ $file -nt $latest ]] && latest=$file
done
echo $latest
python ../../tools/caffemodel2json.py $CAFFE/src/caffe/proto/caffe.proto $latest --data > testNetwork.json
../../tools/json_network_to_rlv_translator.py testNetwork.json > testNetwork.rlv
./computeMNISTAccuracy.py testNetwork.rlv testing_mnist_png
