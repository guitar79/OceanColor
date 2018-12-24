!#/bin/bash

echo "updating aqua LUTS..."
update_luts.py -v aqua


echo "updating terra LUTS..."
update_luts.py -v terra


echo "updating seawifs LUTS..."
update_luts.py -v seawifs


echo "updating viirsn LUTS..."
update_luts.py -v viirsn
