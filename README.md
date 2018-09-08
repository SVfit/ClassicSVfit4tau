# ClassicSVfit4tau
Customized version of SVfit for hh->4tau analysis

# Installation instructions
The SVfitPerformanceStudies package has been tested with CMSSW 9_4_4.
It depends on the following other packages:
	TauAnalysis/ClassicSVfit
	TauAnalysis/SVfitTF

In order to install the code, execute:

```
cmsrel CMSSW_9_4_4
cd CMSSW_9_4_4/src
cmsenv
echo "installing SVfit..."
git clone https://github.com/SVfit/ClassicSVfit4tau TauAnalysis/ClassicSVfit4tau
git clone https://github.com/SVfit/ClassicSVfit TauAnalysis/ClassicSVfit
git clone https://github.com/SVfit/SVfitTF TauAnalysis/SVfitTF
cd $CMSSW_BASE/src
scram b -j 4

echo "installing VAMP library..."
mkdir $CMSSW_BASE/VAMP
wget http://whizard.hepforge.org/oldsrc/vamp-2.3.0.tar.gz -P $CMSSW_BASE/VAMP
tar zxvf $CMSSW_BASE/VAMP/vamp-2.3.0.tar.gz -C $CMSSW_BASE/VAMP
mkdir -p $CMSSW_BASE/VAMP/vamp-2.3.0/prefix/share/doc/vamp
cd $CMSSW_BASE/VAMP/vamp-2.3.0
./configure --prefix=$CMSSW_BASE/VAMP/vamp-2.3.0/prefix
make -j4
make install
cp $CMSSW_BASE/VAMP/vamp-2.3.0/prefix/lib/* $CMSSW_BASE/lib/$SCRAM_ARCH
cp /home/veelken/public/vamp.xml $CMSSW_BASE/config/toolbox/$SCRAM_ARCH/tools/selected/vamp.xml
scram setup vamp
scram tool info vamp
```

In case of compilation problems, please sutmitt an issue on
https://github.com/SVfit/ClassicSVfit4tau/issues

# Running instructions
- [Example(s)](https://github.com/SVfit/ClassicSVfit4tau/blob/master/bin/testClassicSVfit4tau.cc)
