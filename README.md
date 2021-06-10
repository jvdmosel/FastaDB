# FastaDB
Simple database based on [TFClass](http://tfclass.bioinf.med.uni-goettingen.de/)
## Usage
Build image
```console 
docker build -t fastadb-latest . 
```
Run image
```console
docker run --rm -it -v $(pwd)/out:/app/out --name fastadb fastadb-latest
```
