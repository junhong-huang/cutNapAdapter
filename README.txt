# cutNapAdapter

cutNapAdapter: A computational software for cutting adapters of non-capped RNA from NAP-seq data.

Overview:
---------

cutNapAdapter is a software to cut both specific adapters and random barcodes of non-capped RNA from NAP-seq data.

Usage:
---------

Usage:  cutNapAdapter [options] -a <adapter5> -b <adapter3> -1 <read1 fiile, gzip format> -2 <read2 File, gzip format><BR>
[options]<BR>
-o/--output <string>       : output dir[default=cutNapAdapterResults]<BR>
-f/--prefix <string>       : prefix for output files[default=cutNapAdapter]<BR>
-a/--adapter1 <string>     : adapter1 sequence [required]<BR>
-b/--adapter2 <string>     : adapter2 sequence [required for paired-end sequencing]<BR>
-c                         : barcode length of 5'-end[default=6]<BR>
-C                         : barcode length of 3'-end[default=6]<BR>
-m/--min-match             : minimum matches bewteen adapter and sequences[default=5]<BR>
-l/--min-len               : minimum sequence length after trim adapter[default=5]<BR>
-e/--error-rate            : maximum error rate [default<=0.1]<BR>
-v/--verbose               : verbose information<BR>
-V/--version               : cutNapAdapter version<BR>
-h/--help                  : help informations<BR>


Installation:<BR>
---------

Typical install time: within 5 min.<BR>
Download cutNapAdapter-1.0.tar.gz from https://github.com/junhong-huang/cutNapAdapter/releases ; unpack it, and make:<BR>
tar -xzvf cutNapAdapter-1.0.tar.gz<BR>
cd cutNapAdapter-1.0<BR>
make<BR>

System requirements:
---------

Operating system: cutNapAdapter is designed to run on POSIX-compatible platforms, including UNIX, Linux and Mac OS/X. We have tested  most extensively on Linux and MacOS/X because these are the machines we develop on.<BR>
Compiler: The source code is compiled with  the C++ compiler g++. We test the code using the g++ compilers.<BR>
By default, cutNapAdapter does not require any additional libraries to be installed by you.<BR>

run cutNapAdapter:
---------

cutNapAdapter -l 15 -e 0.1 -c 6 -o $outdir -f $prefix -a $adapter1 -b $adapter2 -1 $read1File -2 $read2File<BR>

Output:
---------

The file format of ouput is the same as the input (fq.gz), which can be a direct input as the read mapper STAR or other mappers<BR>

Contact :
---------

*****************************************************************************************<BR>
 \*	cutNapAdapter - A computational software for cutting adapters of non-capped RNA from NAP-seq data.<BR>
 \*<BR>
 \*	Author : Jian-Hua Yang <yangjh7@mail.sysu.edu.cn><BR>
 \* <BR>
 \*	RNA Information Center, School of Life Sciences, Sun Yat-Sen University<BR>
 \*	<BR>
 \*  Create date: 02/07/2022<BR>
 \*  <BR>
 \*  last modified time: 02/07/2022<BR>
 ****************************************************************************************<BR>