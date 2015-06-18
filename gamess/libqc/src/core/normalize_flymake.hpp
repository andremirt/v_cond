
template<class T>
void Basis::normalize(T &A) const {
    for (int b = 0, jj = 0; b < numShells(); ++b) {
	for (int a = 0, ii = 0; a < numShells(); ++a) {

	    for (int j = 0; j < (*this)(b).size(); ++j) {
		for (int i = 0; i < (*this)(a).size(); ++i) {
		    double Ci = (*this)(a).function(i).C;
		    double Cj = (*this)(b).function(j).C;
		    A(ii+i,jj+j) *= Ci*Cj;
		}
	    }
	    ii += (*this)(a).size();
	}
	jj += (*this)(b).size();
    }
}
