#include "spc.hpp"

using namespace std;

SPC::SPC(string &pathname) : Spectra(pathname)
{
    /*************************
    * Reading of an SPC file *
    *************************/
    if(!isOK()) {
        return;
    }
    //cout << "Reading file\n";
	iofile.open(pathname.c_str(), ios_base::in | ios_base::binary);
	if (iofile.is_open()) {
	    //hdr = new SPCHDR;
		// read data as a block:
	    iofile.read((char *)&hdr, sizeof(SPCHDR));
	    // Check the FTFLGS flag
	    if (hdr.ftflgs==TMULTI) {
		    cerr << "Multifile Format\n";
		    exit(2);
	    } else if (hdr.ftflgs==(TMULTI|TORDRD)) {
		    cerr << "Multifile Format Tordrd set\n";
		    exit(2);
	    }
	    // read next data as a block:
        xmin=hdr.ffirst;
	    xmax=hdr.flast;
	    npts=hdr.fnpts;
	    // Read hdr data:
	    shd = new SUBHDR[hdr.fnsub];
	    iofile.read((char *)&shd[0],sizeof(SUBHDR));
	    // Read the Data values (Y axix)
        SetSize(hdr.fnpts);
	    iofile.read((char *)ydata, sizeof(float)*hdr.fnpts);
	    iofile.close();
	}
	// Calculate xdata
	for(int i=0;i<npts;++i) {
        xdata[i] = xmin+i*(xmax-xmin)/(npts-1);
	}
}
SPC::SPC(std::string &filepath, double &hlim, double &llim, int &nps, float *data) : Spectra(filepath, true)
{
    if(!isOK()) {
        return;
    }
    SetSize(nps);
    DWORD one=1;
    //SUBHDR mshd;
    shd = new SUBHDR[one];
    hdr.fnpts=npts;
    hdr.ffirst=hlim;
    hdr.flast=llim;
    hdr.fnsub=one;
    SetYData(data);
    //write(filepath);
    //delete [] shd;
}
void SPC::write()
{
    //string fname = 
    iofile.open(name.c_str(),ios_base::out | ios_base::binary );
    //std::cout << "About to do first write\n";
    iofile.write( (char*)&hdr, sizeof(SPCHDR) );
    //std::cout << "About to do second write\n";
    iofile.write( (char*)&shd[0], sizeof(SUBHDR) );
    //std::cout << "About to do third write\n";
    iofile.write( (char*)ydata, sizeof(float)*hdr.fnpts);
    iofile.close();
    delete [] temp;
}
