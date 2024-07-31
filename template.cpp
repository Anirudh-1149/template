#include <bits/stdc++.h> 
using namespace std;
typedef long long ll; 
typedef vector<int> vi;
typedef pair<int, int> pi;
typedef vector<long long> vll;
#define f(i, a, b) for (int i=a; i<b; i++) 
#define pb push_back
#define mp make_pair

void log(vector<ll> vals)
{
    for(auto val : vals)
        cout<<val<<" ";
    cout<<"\n";
}

ll pow(ll a, ll n, ll m)
{
    ll ans = 1;
    while(n)
    {
        if(n&1)
            ans *=a;
        n = n >>1;
        a = a * a;
        ans%=m;
        a%m;
    }
    return ans;
}


vector<vector<ll>> matrixProd(vector<vector<ll>> &m1, vector<vector<ll>> &m2)
{
    vector<vector<ll>> m3(m1.size(),vector<ll>(m2[0].size(),0));

    f(i,0,m3.size())
    {
        f(j,0,m3[0].size())
        {
            f(k,0,m2.size())
            {
                m3[i][j] += m1[i][k] * m2[k][j];
            }
        }
    }
    return m3;
}

vector<vector<ll>> matrixExponentiation(vector<vector<ll>> &m, ll n)
{
    vector<vector<ll>> ans(m.size(),vector<ll>(m[0].size(),0));
    f(i,0,m.size())
        ans[i][i] = 1;
    while(n)
    {
        if(n&1)
        {
            ans  = matrixProd(ans,m);
        }
        m = matrixProd(m,m);
    }
    return ans;
}

// O(n+m), total number of incremetns of j, cannot exceed the total of increments of i.
vector<int> kmp(string s, string l)
{
    vector<int> ans;
    // ml[i] stores the last index of longest prefix which is also a suffix of string p[0] to p[i]
    vector<int> ml(l.length());
    ml[0] = -1;
    int j = -1;
    f(i,1,l.length())
    {
        while(j != -1 && l[j+1] != l[i])
        {
            j = ml[j];
        }
        if(l[j+1] == l[i])
        {
            j++;
        }
        ml[i] = j;
    }
    j = -1;
    f(i,0,s.length())
    {
        while(j != -1 && s[i] != l[j+1])
        {
            j = ml[j];
        }
        if(s[i] == l[j+1])
        {
            j++;
        }
        if(j == l.length()-1)
        {
            ans.push_back(i);
            j = ml[j];
        }
    }
    return ans;
}

struct Trie
{
    Trie * ch[26];
    bool end;
};

Trie* getNode()
{
    Trie *node = new Trie();
    for(int i = 0;i<26;i++)
    {
        node->ch[i] = nullptr;
    }
    node->end = false;
    return node;
}

Trie * insert(string s, Trie *root)
{
    Trie * temp = root;
    f(i,0,s.length())
    {
        if(temp->ch[s[i]-'a'] == nullptr)
        {
            temp->ch[s[i]-'a'] = getNode();
        }
        temp = temp->ch[s[i]-'a'];
    }
    temp->end = true;
    return temp;
}

bool search(string s, Trie* root)
{
    Trie* temp = root;
    f(i,0,s.length())
    {
        if(temp->ch[s[i]-'a'] == nullptr)   
            return false;
        else 
            temp = temp->ch[s[i]-'a'];
    }
    return temp->end;
}

class SegmentTree{
    private : 
        vector<int> tree;
        int s_;
        int as_;
    void updateValueUtil(int l, int r, int i, int diff, int ti)
    {
        if(i < l || i > r)return;

        tree[ti]+=diff;
        if(l!=r)
        {
            int mid = (l+r)/2;
            updateValueUtil(l,mid,i,diff,ti*2+1);
            updateValueUtil(mid+1,r,i,diff,ti*2+2);
        }
    }
    public :
    SegmentTree(vector<int> &a)
    {
        as_ = a.size();
        s_ = 2*pow(2,(int)ceil(log2(as_))) -1;

        tree.resize(s_);
        constructTree(a,0,as_-1,0);
    }

    int constructTree(vector<int> &a,int l, int r, int ti)
    {
        if(l == r)
        {
            tree[ti] = a[l];
            return a[l];
        }

        int mid = (l+r)/2;
        tree[ti] = constructTree(a,l,mid,ti*2+1) + constructTree(a,mid+1,r,ti*2+2);
        return tree[ti];
    }

    int getSum(int l, int r, int ti, int ql, int qr)
    {
        if(ql <= l && qr >=r)return tree[ti];
        if(r < ql || l > qr) return 0;
        int mid = (l+r)/2;
        return getSum(l,mid,2*ti+1,ql,qr) + getSum(mid+1,r,2*ti+2,ql,qr);
    }

    void updateValue(int i, int diff)
    {
        if(i< 0 || i> as_-1)return;
        updateValueUtil(0,as_-1,i,diff,0);
    }

    
};

class Sat{
    public : 
    int n;
    vector<vector<int>> adj,adjinv;
    vector<bool> used;
    vector<int> order,comp;
    vector<bool> assignment;

    Sat(int n_){
        n = n_;
        adj.resize(2*n);
        adjinv.resize(2*n);
        order.clear();
        used.assign(2*n,false);
        comp.assign(2*n,-1);
    }

    void dfs1(int v)
    {
        used[v] = true;
        for(int u : adj[v])
        {
            if(!used[u])
                dfs1(u);
        }
        order.push_back(v);
    }

    void dfs2(int v, int c)
    {
        comp[v] = c;
        for(int u : adjinv[v])
        {
            if(comp[u] == -1)
                dfs2(u,c);
        }
    }

    bool solve()
    {
        for(int i = 0;i<2*n;i++)
            if(!used[i])
                dfs1(i);
        
        for(int i = 0,j = 0;i<n;i++)
        {
            int v = order[2*n-i-1];
            if(comp[v] == -1)
            {
                dfs2(v,j++);
            }
        }

        assignment.assign(n,false);

        for(int i = 0;i<2*n;i+=2)
        {
            if(comp[i]==comp[i+1])
                return false;
            assignment[i/2] = comp[i] > comp[i+1];
        }
        return true;
    }

    void add_disjunction(int a, bool na, int b, bool nb)
    {
        a = 2*a ^ na;
        b = 2*b ^ nb;
        int nega = a^1;
        int negb = b^1;
        adj[nega].push_back(b);
        adj[negb].push_back(a);
        adjinv[a].push_back(negb);
        adjinv[b].push_back(nega);
    }
};


vector<int> manachers(string s)
{
    int n = s.length();
    string st = "$";
    for(char c: s)
    {
        st.push_back(c);
        st.push_back('$');
    }

    vector<int> lps(st.length(),0);
    int start = 0;
    int end = 0;
    for(int i = 0;i<st.length();)
    {
        while(start > 0 && end < st.length()-1 && st[start-1] == st[end+1])
        {
            start--;
            end++;
        }

        lps[i] = end-start+1;
        
        int nc = end + 1;

        for(int j = i+1;j<=end;j++)
        {
            lps[j] = min(lps[2*i-j], 2*(end-j)+1);
            if(j + lps[j]/2 == end && end != st.length()-1)
            {
                nc = j;
                break;
            }
        }

        if(end == st.length()-1)break;
        i = nc;
        end = i + lps[i]/2;
        start = i - lps[i]/2;
    }
    return lps;
}

class Lca{
    public : 
    int n,l;
    vector<vector<int>> adj;
    int timer;
    vector<int> tin, tout;
    vector<vector<int>> up;

    void dfs(int v, int p)
    {
        tin[v] = ++timer;
        up[v][0] = p;
        for(int i = 1;i<=l;i++)
            up[v][i] = up[up[v][i-1]][i-1];
        
        for(int u : adj[v])
        if(u != p)
            dfs(u,v);
        
        tout[v] = ++timer;
    }

    bool is_ancestor(int u, int v)
    {
        return tin[u] <= tin[v] && tout[u] >= tout[v];
    }

    int lca(int u, int v)
    {
        if(is_ancestor(u,v))return u;
        if(is_ancestor(v,u)) return v;

        for(int i = l;i>=0;i--)
        {
            if(!is_ancestor(up[u][i],v))
                u = up[u][i];
        }
        return up[u][0];
    }

    void preprocess(int root)
    {
        tin.resize(n);
        tout.resize(n);
        timer = 0;
        l = ceil(log2(n));
        up.assign(n,vector<int>(l+1));
        dfs(root,root);
    }
};

// clockwise orientation means that area is less than 0.
int orientation(vector<int> a, vector<int> b, vector<int> c)
{
    double v = a[0] * (b[1] - c[1]) + b[0] * (c[1] - a[1]) + c[0] * (a[1] - b[1]);
    if(v < 0) return -1;
    if(v > 0)return 1;
    return 0;
}

// return if the points are clockwise
bool cw(vector<int> a, vector<int> b, vector<int> c)
{
    int o = orientation(a,b,c);
    return o < 0;
}

vector<vector<int>> convexHull(vector<vector<int>> &points)
{
    vector<int> p0 = *min_element(points.begin(),points.end(), [](vector<int> a, vector<int> b){
        return a < b;
    });

    // sorting in clockwise direction
    sort(points.begin(),points.end(),[&p0](const vector<int> &a, vector<int> &b){
        int o = orientation(p0,a,b);
        if(o == 0)
            return pow(p0[0]-a[0],2) + pow(p0[1] - a[1],2) < pow(p0[0] - b[0],2) + pow(p0[1] - b[1],2);
        return o < 0;
    });

    vector<vector<int>> st;
    st.push_back(points[0]);
    st.push_back(points[1]);

    for(int i = 2;i<points.size();i++)
    {
        if(st.size() > 1 && !cw(st[st.size()-2],st[st.size()-1],points[i]))
            st.pop_back();
        st.push_back(points[i]);
    }

    return st;
}

// z[i] is the length of the string that is prefix of string s, also a prefix of the string s[i,n-1]
vector<int> z_function(string s) {
    int n = s.size();
    vector<int> z(n);
    int l = 0, r = 0;
    for(int i = 1; i < n; i++) {
        if(i < r) {
            z[i] = min(r - i, z[i - l]);
        }
        while(i + z[i] < n && s[z[i]] == s[i + z[i]]) {
            z[i]++;
        }
        if(i + z[i] > r) {
            l = i;
            r = i + z[i];
        }
    }
    return z;
}

class Bridge{
    public :
    vector<vector<int>> adj;
    Bridge(vector<vector<int>> &adj_)
    {
        adj = adj_;
    }

    vector<pair<int,int>> find_bridges()
    {
        int n = adj.size();
        vector<int> tin(n);
        vector<bool> visi(n);
        vector<int> low(n);
        vector<pair<int,int>> ans;
        dfs(0,-1,tin,low,visi,0,ans);
        return ans;
    }

    void dfs(int node, int p, vector<int> &tin,vector<bool>&visi,vector<int>&low, int &timer,vector<pair<int,int>> &ans)
    {
        visi[node] = true;
        tin[node] = low[node] = timer++;
        bool skip_parent = false;

        for(int v : adj[node])
        {
            if(v == p && !skip_parent)
            {
                skip_parent = true;
                continue;
            }
            if(visi[v])
            {
                low[node] = min(low[node],tin[v]);
            }
            else 
            {
                dfs(v,node,tin,visi,low,timer,ans);
                low[node] = min(low[node],low[v]);
                if(low[v] > tin[node])
                {
                    ans.push_back({v,node});
                }
            }
        }
    }
};

void solve()
{
    vector<vector<int>> points = {{0,3},{0,0},{3,0},{3,3},{2,2}};
    for(auto p : convexHull(points))
    {
        cout<<p[0]<< " "<<p[1]<<"\n";
    }
}


int main()
{
    ios::sync_with_stdio(false); cin.tie(0); cout.tie(0);
    int t = 1;
    // cin>>t;
    while(t--)
    {
        solve();
    }
    return 0;
}
