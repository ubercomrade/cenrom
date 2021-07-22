# cenrom
Calculate enrichment of matrix in data

## Requirements

PYTHON:
  * numpy: `pip3 install numpy`
  * MOODS (https://www.cs.helsinki.fi/group/pssmfind/,
  https://github.com/jhkorhonen/MOODS): `pip3 install MOODS-python`

## Installation

```  
git clone https://github.com/ubercomrade/pipeline.git  
cd pipeline/  
pip3 install -e .  
```

## Usage
The command `pipeline.py -h` return:

```
usage: cenrom.py [-h] [-n NTIMES] [-f N] [-t THRESHOLD] fasta matrix N

positional arguments:
  fasta                 path to FASTA file
  matrix                path to MATRIX file (pcm or pfm depend on source)
  N                     promoters of organism (hg38, mm10, tair10)

optional arguments:
  -h, --help            show this help message and exit
  -n NTIMES, --ntimes NTIMES
                        N times FASTA for background, def=100
  -f N, --format N      [homer, cisbp, hocomoco] format of input PFM (HOMER,
                        CISBP) or PCM (HOCOMOCO)
  -t THRESHOLD, --threshold THRESHOLD
                        threshold based on FPR, def=1.9*10^(-4)
```
### Example run
```
cenrom.py peaks.fa ARNT2_MOUSE.H11MO.0.D.pcm mm10 -n 50 -f hocomoco 
```

### Positional arguments description

**First positional argument**:
```
fasta                   path to fasta file
```
You have to give path to fasta file. 

**Second positional argument**:
```
matrix                     path to MATRIX file
```
You have to give path to matrix file. The format of file deneond on `-f N, --format N` argument

**Third positional argument**:
```
N                promoters of organism (hg38, mm10, tair10)
```
Value of N can be _hg38_ or _mm10_ or _tair10_. Depend on organism used in research

### Flag arguments description

**First flag**:

```
-n NTIMES, --ntimes NTIMES                N times FASTA for background, def=100
```
How many times the background will be larger than the original FASTA file

**Second flag**:


```
-f N, --format N      [homer, cisbp, hocomoco]
```
Choose format of your matrix.

### Format examples
**homer**:
```
>RTWCAATGWATC	1-RTWCAATGWATC,BestGuess:EIL4(EIL)/Tomato-EIL4-ChIP-Seq(GSE116581)/Homer(0.983)	8.109841	-1085.337425	0	T:449.0(55.85%),B:1272.8(2.69%),P:1e-471
0.538	0.001	0.460	0.001
0.001	0.001	0.001	0.997
0.386	0.052	0.110	0.451
0.066	0.932	0.001	0.001
0.848	0.001	0.150	0.001
0.471	0.104	0.158	0.266
0.001	0.126	0.001	0.872
0.001	0.001	0.983	0.015
0.486	0.137	0.034	0.344
0.997	0.001	0.001	0.001
0.001	0.343	0.001	0.655
0.160	0.441	0.150	0.250
```

**cisbp**:
```
Pos	A	C	G	T
1	0.328843577877036	0.158472676531907	0.193558964202559	0.319124781388498
2	0.340062470089621	0.138259082775898	0.138259082775898	0.383419364358583
3	0.250780641145707	0.152106147730106	0.248285332074236	0.348827879049951
4	0.146662255535782	0.212328023990482	0.47729386085815	0.163715859615587
5	0.0491484979435091	0.10938199373243	0.841214166900551	0.00025534142351
6	0.0986661116261265	0.881602734252051	9.02047344252513e-06	0.0197221336483805
7	0.0012280856392551	1.33726536030681e-05	0.977948064627909	0.0208104770792331
8	0.0225806909007298	0.624589825817514	0.326064480811723	0.0267650024700323
9	0.13556626962152	0.284138159031754	0.361210393228366	0.21908517811836
10	0.23482140428296	0.292605487709977	0.182879093545217	0.289694014461846
```

**HOCOMOCO**
```
>ATF4_HUMAN.H11MO.0.A
199.00000000000003	32.00000000000001	191.00000000000003	78.00000000000001
118.00000000000001	93.00000000000001	197.00000000000003	92.00000000000001
267.00000000000006	165.00000000000003	68.00000000000001	0.0
1.0000000000000002	1.0000000000000002	0.0	498.00000000000006
0.0	3.000000000000001	489.00000000000006	8.000000000000002
490.00000000000006	3.000000000000001	5.000000000000001	2.0000000000000004
0.0	26.000000000000004	3.000000000000001	471.00000000000006
0.0	3.000000000000001	497.00000000000006	0.0
51.00000000000001	401.00000000000006	1.0000000000000002	47.00000000000001
497.00000000000006	3.000000000000001	0.0	0.0
499.00000000000006	0.0	1.0000000000000002	0.0
8.000000000000002	146.00000000000003	65.00000000000001	281.00000000000006
```

**Third flag**:

```
-t THRESHOLD, --threshold THRESHOLD
```
The argument `-t/--threshold` sets up FPR for choosing threshold score based on FPR value [1]. The default value is equal to 0.00019.

## Useful links

 * [cisbp](http://cisbp.ccbr.utoronto.ca/)
 * [hocomoco](https://hocomoco11.autosome.ru/)

## Reference

 [1] Victor Levitsky, Elena Zemlyanskaya, Dmitry Oshchepkov, Olga Podkolodnaya, Elena Ignatieva, Ivo Grosse, Victoria Mironova, Tatyana Merkulova, A single ChIP-seq dataset is sufficient for comprehensive analysis of motifs co-occurrence with MCOT package, Nucleic Acids Research, Volume 47, Issue 21, 02 December 2019, Page e139, https://doi.org/10.1093/nar/gkz800

