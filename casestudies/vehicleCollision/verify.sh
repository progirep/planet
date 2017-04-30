set -e
unset -v latest
for file in "snapshots"/_iter_*.caffemodel; do
  [[ $file -nt $latest ]] && latest=$file
done
echo $latest
python ../../tools/caffemodel2json.py $CAFFE/src/caffe/proto/caffe.proto $latest --data > testNetworkWithPooling.json
./json_net_against_csv_robustness_verifier.py  testNetworkWithPooling.json collisions.csv
