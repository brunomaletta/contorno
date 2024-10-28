#include <bits/stdc++.h>

using namespace std;

#define _ ios_base::sync_with_stdio(0);cin.tie(0);
#define endl '\n'

typedef long long ll;

const int INF = 0x3f3f3f3f;
const ll LINF = 0x3f3f3f3f3f3f3f3fll;

using namespace chrono;
struct timer : high_resolution_clock {
	const time_point start;
	timer(): start(now()) {}
	int operator()() {
		return duration_cast<milliseconds>(now() - start).count();
	}
};

typedef double ld;
const ld DINF = 1e18;
const ld pi = acos(-1.0);
const ld eps = 1e-9;

#define sq(x) ((x)*(x))

bool eq(ld a, ld b) {
	return abs(a - b) <= eps;
}

struct pt { // ponto
	ld x, y;
	pt(ld x_ = 0, ld y_ = 0) : x(x_), y(y_) {}
	bool operator < (const pt p) const {
		if (!eq(x, p.x)) return x < p.x;
		if (!eq(y, p.y)) return y < p.y;
		return 0;
	}
	bool operator == (const pt p) const {
		return eq(x, p.x) and eq(y, p.y);
	}
	pt operator + (const pt p) const { return pt(x+p.x, y+p.y); }
	pt operator - (const pt p) const { return pt(x-p.x, y-p.y); }
	pt operator * (const ld c) const { return pt(x*c  , y*c  ); }
	pt operator / (const ld c) const { return pt(x/c  , y/c  ); }
	ld operator * (const pt p) const { return x*p.x + y*p.y; }
	ld operator ^ (const pt p) const { return x*p.y - y*p.x; }
	friend istream& operator >> (istream& in, pt& p) {
		return in >> p.x >> p.y;
	}
};

ld angle(pt v) { // angulo do vetor com o eixo x
	ld ang = atan2(v.y, v.x);
	if (ang < 0) ang += 2*pi;
	return ang;
}

namespace blossom{
	#define d(x) (lab[x.u]+lab[x.v]-e[x.u][x.v].w*2)
	const int N=403*2; using ll = long long; using T = int; // sum of weight, single weight
	const T inf=numeric_limits<T>::max()>>1;
	struct Q{ int u, v; T w; } e[N][N]; vector<int> p[N];
	int n, m=0, id, h, t, lk[N], sl[N], st[N], f[N], b[N][N], s[N], ed[N], q[N]; T lab[N];
	void upd(int u, int v){ if(!sl[v] || d(e[u][v]) < d(e[sl[v]][v])) sl[v] = u; }
	void ss(int v){
		sl[v]=0; for(auto u=1; u<=n; u ++) if(e[u][v].w > 0 && st[u] != v && !s[st[u]]) upd(u, v);
	}
	void ins(int u){ if(u <= n) q[++ t] = u; else for(auto v : p[u]) ins(v); }
	void mdf(int u, int w){ st[u]=w; if(u > n) for(auto v : p[u]) mdf(v, w); }
	int gr(int u,int v){
		if((v=find(p[u].begin(), p[u].end(), v) - p[u].begin()) & 1){
			reverse(p[u].begin()+1, p[u].end()); return (int)p[u].size() - v;
		}
		return v;
	}
	void stm(int u, int v){
		lk[u] = e[u][v].v;
		if(u <= n) return; Q w = e[u][v];
		int x = b[u][w.u], y = gr(u,x);
		for(auto i=0; i<y; i ++) stm(p[u][i], p[u][i^1]);
		stm(x, v); rotate(p[u].begin(), p[u].begin()+y, p[u].end());
	}
	void aug(int u, int v){
		int w = st[lk[u]]; stm(u, v); if(!w) return;
		stm(w, st[f[w]]); aug(st[f[w]], w);
	}
	int lca(int u, int v){
		for(++ id; u|v; swap(u, v)){
			if(!u) continue; if(ed[u] == id) return u;
			ed[u] = id; if(u = st[lk[u]]) u = st[f[u]]; // not ==
		}
		return 0;
	}
	void add(int u, int a, int v){
		int x = n+1; while(x <= m && st[x]) x ++;
		if(x > m) m ++;
		lab[x] = s[x] = st[x] = 0; lk[x] = lk[a];
		p[x].clear(); p[x].push_back(a);
		for(auto i=u, j=0; i!=a; i=st[f[j]]) p[x].push_back(i), p[x].push_back(j=st[lk[i]]), ins(j);
		reverse(p[x].begin()+1, p[x].end());
		for(auto i=v, j=0; i!=a; i=st[f[j]]) p[x].push_back(i), p[x].push_back(j=st[lk[i]]), ins(j);
		mdf(x, x); for(auto i=1; i<=m; i ++) e[x][i].w = e[i][x].w = 0;
		memset(b[x]+1, 0, n*sizeof b[0][0]);
		for(auto u : p[x]){
			for(v=1; v<=m; v ++) if(!e[x][v].w || d(e[u][v]) < d(e[x][v])) e[x][v] = e[u][v],e[v][x] = e[v][u];
			for(v=1; v<=n; v ++) if(b[u][v]) b[x][v] = u;
		}
		ss(x);
	}
	void ex(int u){  // s[u] == 1
		for(auto x : p[u]) mdf(x, x);
		int a = b[u][e[u][f[u]].u],r = gr(u, a);
		for(auto i=0; i<r; i+=2){
			int x = p[u][i], y = p[u][i+1];
			f[x] = e[y][x].u; s[x] = 1; s[y] = 0; sl[x] = 0; ss(y); ins(y);
		}
		s[a] = 1; f[a] = f[u];
		for(auto i=r+1; i<p[u].size(); i ++) s[p[u][i]] = -1, ss(p[u][i]);
		st[u] = 0;
	}
	bool on(const Q &e){
		int u=st[e.u], v=st[e.v], a;
		if(s[v] == -1) f[v] = e.u, s[v] = 1, a = st[lk[v]], sl[v] = sl[a] = s[a] = 0, ins(a);
		else if(!s[v]){
			a = lca(u, v); if(!a) return aug(u,v), aug(v,u), true; else add(u,a,v);
		}
		return false;
	}
	bool bfs(){
		memset(s+1, -1, m*sizeof s[0]); memset(sl+1, 0, m*sizeof sl[0]);
		h = 1; t = 0; for(auto i=1; i<=m; i ++) if(st[i] == i && !lk[i]) f[i] = s[i] = 0, ins(i);
		if(h > t) return 0;
		while(true){
			while(h <= t){
				int u = q[h ++];
				if(s[st[u]] != 1) for(auto v=1; v<=n; v ++) if(e[u][v].w > 0 && st[u] != st[v])
					if(d(e[u][v])) upd(u, st[v]); else if(on(e[u][v])) return true;
			}
			T x = inf;
			for(auto i=n+1; i<=m; i ++) if(st[i] == i && s[i] == 1) x = min(x, lab[i]>>1);
			for(auto i=1; i<=m; i ++) if(st[i] == i && sl[i] && s[i] != 1) x = min(x, d(e[sl[i]][i])>>s[i]+1);
			for(auto i=1; i<=n; i ++) if(~s[st[i]]) if((lab[i] += (s[st[i]]*2-1)*x) <= 0) return false;
			for(auto i=n+1 ;i<=m; i ++) if(st[i] == i && ~s[st[i]]) lab[i] += (2-s[st[i]]*4)*x;
			h = 1; t = 0;
			for(auto i=1; i<=m; i ++) if(st[i] == i && sl[i] && st[sl[i]] != i && !d(e[sl[i]][i]) && on(e[sl[i]][i])) return true;
			for(auto i=n+1; i<=m; i ++) if(st[i] == i && s[i] == 1 && !lab[i]) ex(i);
		}
		return 0;
	}
	template<typename TT> pair<ll, vector<array<int, 2>>> run(int N, vector<tuple<int,int,TT>> edges){ // 1-based
		for(auto &[u, v, w]: edges) ++ u, ++ v;
		memset(ed+1, 0, m*sizeof ed[0]); memset(lk+1, 0, m*sizeof lk[0]);
		n = m = N; id = 0; iota(st+1, st+n+1, 1); T wm = 0; ll weight = 0;
		for(auto i=1; i<=n; i ++) for(auto j=1; j<=n; j ++) e[i][j] = {i,j,0};
		for(auto [u,v,w] : edges) wm = max(wm, e[v][u].w=e[u][v].w=max(e[u][v].w,(T)w));
		for(auto i=1; i<=n; i ++) p[i].clear();
		for(auto i=1; i<=n; i ++) for(auto j=1; j<=n; j ++) b[i][j] = i*(i==j);
		fill_n(lab+1, n, wm); while(bfs());
		vector<array<int, 2>> matching;
		for(auto i=1; i<=n; i ++) if(i < lk[i]) weight += e[i][lk[i]].w, matching.push_back({i - 1, lk[i] - 1});
		return {weight, matching};
	}
	#undef d
}

const int MAX = 5000;

ll d[MAX][MAX];
int prv[MAX][MAX];
int n;
vector<pt> coord;
map<int, int> mp, orig;

void floyd_warshall() {
	for (int k = 0; k < n; k++) for (int i = 0; i < n; i++) for (int j = 0; j < n; j++)
		if (d[i][k] + d[k][j] < d[i][j]) {
			d[i][j] = d[i][k] + d[k][j];
			prv[i][j] = prv[k][j];
		}
}

vector<pair<int, int>> g[MAX]; // v, id
int used[MAX];
int vis[MAX];

int main() { _
	string line;
	bool started = false;
	vector<tuple<int, int, int>> edg;
	while (getline(cin, line)) {
		if (!line.size()) {
			started = true;
			continue;
		}
		if (!started) {
			stringstream ss(line);
			double a, b;
			ss >> a >> b;
			coord.emplace_back(a, b);
			continue;
		}
		int a, b;
		double c;
		stringstream ss(line);
		ss >> a >> b >> c;
		edg.emplace_back(a, b, (int) c);
	}
	auto get_id = [&](int i) {
		auto it = mp.find(i);
		if (it != mp.end()) return it->second;
		int sz = mp.size();
		mp[i] = sz;
		orig[sz] = i;
		return sz;
	};
	for (int i = 0; i < MAX; i++) {
		memset(d[i], INF, sizeof(d[i]));
		memset(prv[i], -1, sizeof(prv[i]));
	}
	for (int i = 0; i < MAX; i++) {
		d[i][i] = 0;
		prv[i][i] = i;
	}
	set<pair<int, int>> st;
	vector<tuple<int, int, int>> edg2;
	for (auto& [a, b, c] : edg) {
		a = get_id(a);
		b = get_id(b);
		if (a > b) swap(a, b);
	}
	for (auto& [a, b, c] : edg) {
		if (st.count(pair<int, int>(a, b))) continue;
		st.emplace(a, b);
		edg2.emplace_back(a, b, c);
	}
	edg = edg2;
	ll sum = 0;
	for (auto& [a, b, c] : edg) {
		sum += c;
		d[a][b] = d[b][a] = c;
		prv[a][b] = a;
		prv[b][a] = b;
		n = max({n, a, b});
	}
	n++;
	cerr << "n: " << n << endl;
	cerr << "m: " << edg.size() << endl;;
	cerr << "sum of edges: " << sum << endl;
	vector<int> deg(n);
	for (auto [a, b, c] : edg) deg[a]++, deg[b]++;

	floyd_warshall();
	vector<tuple<int, int, ll>> blossom_edge;
	for (int i = 0; i < n; i++) if (deg[i] % 2)
		for (int j = i+1; j < n; j++) if (deg[j] % 2)
			blossom_edge.emplace_back(i, j, 1000000-d[i][j]);
	auto [weight, matching] = blossom::run(n, blossom_edge);
	weight = 1000000*matching.size() - weight;

	mt19937 rng(0);
	//shuffle(edg.begin(), edg.end(), rng);

	ll ans = sum;
	for (auto [a, b] : matching) {
		while (a != b) {
			int a2 = prv[b][a];
			deg[a]++, deg[a2]++;
			edg.emplace_back(a, a2, d[a][a2]);
			ans += d[a][a2];
			assert(a2 != -1);
			a = a2;
		}
	}

	cerr << "final ans: " << ans << endl;

	for (int i = 0; i < n; i++) {
		assert(deg[i]%2 == 0);
		assert(deg[i] > 0);
	}

	for (int i = 0; i < edg.size(); i++) {
		auto [a, b, c] = edg[i];
		g[a].emplace_back(b, i);
		g[b].emplace_back(a, i);
	}

	int INI = 0;
	vector<int> tour = {INI};
	int at = INI;
	bool lap = false;
	int num_corners = 0, num_straights = 0;
	while (tour.size() < edg.size()+1) {
		auto ang = [&](int i) -> double {
			if (tour.size() < 2) return 0;
			double v1 = angle(coord[orig[i]] - coord[orig[tour.end()[-1]]]);
			double v2 = angle(coord[orig[tour.end()[-1]]] - coord[orig[tour.end()[-2]]]);
			double th = abs(v1 - v2);
			if (th > pi) th = 2*pi - th;
			if (!lap) {
				th = v1 - v2 + pi/2 + 2*pi;
				while (th > 2*pi) th -= 2*pi;
			}
			return th;
		};
		sort(g[at].begin(), g[at].end(), [&](auto a, auto b) {
			double th1 = ang(a.first), th2 = ang(b.first);
			if (abs(th1 - th2) < eps) return false;
			return th1 < th2;
		});
		for (auto [i, id] : g[at]) if (!used[id]) {
			deg[at]--, deg[i]--;
			used[id] = 1;
			for (int j = 0; j < n; j++) vis[j] = 0;

			int fi = 0;
			while (deg[fi] == 0) fi++;
			vector<int> st = {fi};
			vis[fi] = 1;
			while (st.size()) {
				int j = st.back(); st.pop_back();
				for (auto [k, id2] : g[j]) if (!used[id2] and !vis[k]) {
					vis[k] = 1;
					st.push_back(k);
				}
			}
			bool con = true;
			for (int j = 0; j < n; j++) if (!vis[j] and deg[j] > 0) con = false;
			int cnt_odd = 0;
			for (int j = 0; j < n; j++) if (deg[j]%2) cnt_odd++;
			if (con and cnt_odd <= 2) {
				at = i;
				tour.push_back(at);
				if (at == INI) lap = true;

				double v1 = angle(coord[orig[i]] - coord[orig[tour.end()[-1]]]);
				double v2 = angle(coord[orig[tour.end()[-1]]] - coord[orig[tour.end()[-2]]]);
				double th = abs(v1 - v2);
				if (th > pi) th = 2*pi - th;
				if (lap and th > pi/12) num_corners++;
				else num_straights++;
				break;
			}
			used[id] = 0;
			deg[at]++, deg[i]++;
		}
	}

	assert(tour[0] == tour.end()[-1]);

	reverse(tour.begin(), tour.end());

	cerr << "number of corners: " << num_corners << endl;
	cerr << "number of straights: " << num_straights << endl;

	cout << fixed << setprecision(8);
	cout << n << endl;
	for (int i = 0; i < n; i++) cout << coord[orig[i]].x << " "
		<< coord[orig[i]].y << endl;
	for (auto i : tour) cout << i << " ";
	cout << endl;
	exit(0);
}
