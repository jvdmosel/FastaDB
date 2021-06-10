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
## Examples:
Output non-aligned sequences for given species:
``` console 
Gorilla gorilla gorilla
```
``` console 
Homo sapiens
```
Output specific TFClass node:
``` console 
1.2.2.2
```
Output aligned TFClass node:
``` console 
1.2 -a
```
Compare the shannon-entropy of an aligned level 4 node with a level 5 node:
``` console 
-comp 2.1.3.1 2.1.3.1.1
```
