/**
 * @file src/edge-detect/nelder_mead.inl
 * @author Thomas Verduin <T.Verduin@tudelft.nl>
 * @author Sebastiaan Lokhorst <S.R.Lokhorst@tudelft.nl>
 */

#ifndef CCPO__NELDER_MEAD__INLINE_INCLUDED
#define CCPO__NELDER_MEAD__INLINE_INCLUDED

namespace nelder_mead {

template<int N> simplex<N> randomize(const params<N>& min, const params<N>& max) {
	simplex<N> simplex;
	for(int i = 0; i < N+1; i++)
		for(int j = 0; j < N; j++)
			simplex[i][j] = min[j] + (max[j]-min[j])*drand48();
	return simplex;
}

template<int N> void bound(simplex<N>& simplex, const params<N>& min, const params<N>& max) {
	for(int i = 0; i < N+1; i++)
		for(int j = 0; j < N; j++) {
			if(simplex[i][j] <= min[j]) {
				simplex = randomize<N>(min, max);
				return;
			} else if(simplex[i][j] >= max[j]) {
				simplex = randomize<N>(min, max);
				return;
			}
		}
}

template<int N> double iterate(const function<N>& fun, simplex<N>& simplex) {
	std::array<double, N+1> feval;
	for(int i = 0; i < N+1; i++)
		feval[i] = fun(simplex[i]);
	for(int j = 1; j < N+1; j++) {
		const auto f = feval[j];
		const auto v = simplex[j];
		int i = j-1;
		while((i > -1) && (feval[i] > f)) {
			feval[i+1] = feval[i];
			simplex[i+1] = simplex[i];
			i--;
		}
		feval[i+1] = f;
		simplex[i+1] = v;
	}
	params<N> vertex_g = {0};
	for(int i = 0; i < N; i++)
		for(int j = 0; j < N; j++)
			vertex_g[j] += simplex[i][j]/N;
	params<N> vertex_r;
	params<N> vertex_e;
	params<N> vertex_c;
	for(int j = 0; j < N; j++) {
		vertex_r[j] = vertex_g[j] + 1.0d*(vertex_g[j]-simplex[N][j]);
		vertex_e[j] = vertex_g[j] + 2.0d*(vertex_g[j]-simplex[N][j]);
		vertex_c[j] = vertex_g[j] - 0.5d*(vertex_g[j]-simplex[N][j]);
	}
	const double feval_r = fun(vertex_r);
	const double feval_e = fun(vertex_e);
	const double feval_c = fun(vertex_c);
	if((feval_r >= feval[0]) && (feval_r < feval[N-1])) {
		simplex[N] = vertex_r;
	} else if(feval_r < feval[0]) {
		if(feval_e < feval_r)
			simplex[N] = vertex_e;
		else
			simplex[N] = vertex_r;
	} else if(feval_c < feval[0]) {
		simplex[N] = vertex_c;
	} else {
		for(int i = 1; i < N+1; i++)
			for(int j = 0; j < N; j++)
				simplex[i][j] = simplex[0][j] + 0.5d*(simplex[i][j]-simplex[0][j]); 
	}
	double eps = 0;
	for(int j = 0; j < N; j++) {
		const double err = fabs(simplex[N][j]-simplex[0][j])/fabs(simplex[0][j]);
		eps = std::max(eps, err);
	}
	return eps;
}

}

#endif