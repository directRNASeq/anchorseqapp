# Direct Sequencing of tRNA-Phe using Anchor-based Algorithm

Copyright (c) 2020 Wenjia Li, Shenglong Zhang, New York Institute of Technology

This software is licensed under the the MIT License. Contact for other licensing 
terms.

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.

## General Introduction

This software provides a proof-of-concept implementation of an anchor-based algorithm 
that can de novo and directly sequence RNA samples using MFE data output from LC-MS.


### Requirements

To successfully install and run this software, you should download and install the
following software.

1. Java SE Development Kit (JDK) version "12.0.1" or above
2. Maven
3. Eclipse Jetty Server (https://www.eclipse.org/jetty/download.html)
4. a Java IDE (Either Eclipse or IntelliJ)


### Installation

1. Download and unzip Eclipse Jetty Server. You should see a directory named 
`jetty-distribution-9.4.28.v20200408`.

2. Download and install the appropriate JDK, go to the root directory of the project,
and run the following command in the terminal: 

```bash
mvn package
cp target/sequencing-0.0.1-SNAPSHOT.war jetty-distribution-9.4.28.v20200408/webapps/
cp seq_web_app/sequencing_app.html jetty-distribution-9.4.28.v20200408/webapps/
cp seq_web_app/plot_mass.py  jetty-distribution-9.4.28.v20200408/webapps/
mkdir jetty-distribution-9.4.28.v20200408/config/
cp config/* jetty-distribution-9.4.28.v20200408/config/
```


### Usage of This Software in Java IDE

To run the software directly from Java, you can open the project from an IDE.

The FindSequence.java file in the software contains a main method where you can execute it. In the
main method, there is a statement as follows.

`String fileName = "/data/input_data/TableS1_052118s04.txt";`

You can replace the fileName with any input data file that has the same format as the sample data files 
(which are stored in the `data/input_data` directory). Note that the full file path is required to allow 
the software locate the data file and then read out sequence successfully.

In the IDE such as Eclipse and IntelliJ, you can simply navigate to the FindSequence.java file, and
click the "Run" button or menu item to run the whole software.

### Usage of This Software in Web Interface

1. To run the software from web interface, you should first ensure that the `sequencing-0.0.1-SNAPSHOT.war`
file is located at `jetty-distribution-9.4.28.v20200408/webapps` directory. Then run the following commands.

```bash
cd jetty-distribution-9.4.28.v20200408
java -jar start.jar
```

2. Set the anchor with the `anchor_bank.csv` which is located at `jetty-distribution-9.4.28.v20200408/config` directory. Please specify 
only one anchor at each time, in order to achieve the best reading accuracy. 

3. Open web-based sequencing app which is located at `jetty-distribution-9.4.28.v20200408/webapps/sequencing_app.html`.
Select and upload the data file which is located at `data/input_data`. 

4. The results will appear in the `jetty-distribution-9.4.28.v20200408/Result` and `jetty-distribution-9.4.28.v20200408/fastaResult` 
directories. The result can also be visualized as a RT-Mass plot in the web interface.


### Parameters

`<anchor_bank.csv>`: The anchor dataset filename which is located at `config` directory
	Filename of a list of anchor in the format of `AnchorName<TabSpace>AnchorMass`

`<base_bank.csv>`: The base and modification dataset filename which is located at `config` directory
	Filename of a list of bases and their known modifications in the format of `BaseName<TabSpace>BaseMass`

`<directory>`: The directory to output the result for final sequence reads

`<fasta_directory>`: The directory to output the result in terms of the FASTA format

### Input Data

The input data (in the txt format and located at the `data/input_data` directory) is generated based 
on the MFE data file (in the xls format) that is produced by the Agilent MassHunter Qualitative 
Analysis software. In the input data file, you should ensure that it is a list of compounds with 
one compound per row and contains the following columns for each compound:

`Mass<TabSpace>RT<TabSpace>Vol<TabSpace>Cpd<TabSpace>QS`

These columns contain the following information:

Mass: neutral mass RT: retention time Vol: integrated intensity Cpd: compound ID (generally an integer) QS: Quality Score

Processing will fail if any of these columns are missing, or if they are in the wrong order.


