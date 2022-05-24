#include <fstream>
#include <iostream>
#include <string>
#include <vector>
#include <sstream>
#include <algorithm>

//klasa z sekwencja
class SEQ {
public:
    std::vector<int> qual;
    std::vector<char> fasta;
    std::vector<int> pos;
};

//pobieranie 5 roznych sekwencji z pliku
std::vector<SEQ> getInstance(int instanceNum){
    std::vector<SEQ> instance(5);
    std::ifstream ffasta("instancje/" + std::to_string(instanceNum) + "/inst.fasta");
    std::ifstream fqual("instancje/"+ std::to_string(instanceNum) + "/inst.qual");

    char nuc;
    int seq_num = -1;
    int nuc_num = 1;
    std::string ignore;
//czytanie fasta
    if(ffasta.is_open()) {
        do{
            nuc = ffasta.get();
            if(nuc == EOF) break;
            if (nuc == '>'){
                std::getline(ffasta, ignore);
                seq_num++;
                nuc_num = 1;
                continue;
            } else if(nuc == '\n'){
                continue;
            } else {
                instance[seq_num].fasta.push_back(nuc);
                instance[seq_num].pos.push_back(nuc_num);
                nuc_num++;
            }
        } while (nuc != EOF);
    } else {
        std::cout<< "Nie znaleziono pliku FASTA!";
        std::cin.ignore();
        
        exit(1);
    }

//czytanie qual
    std::string q;
    seq_num = -1;
    if(fqual.is_open()) {
        while (fqual >> q) {
            if (q[0] == '>') {
                seq_num++;
                std::getline(ffasta, ignore);
                continue;
            } else if (q == "\n") {
                continue;
            } else {
                std::stringstream toNum(q);
                int num;
                toNum >> num;
                if (num == 0) continue;
                instance[seq_num].qual.push_back(num);
            }
        }
    } else {
        std::cout<< "Nie znaleziono pliku QUAL!";
        std::cin.ignore();
        exit(2);
    }

    return instance;
}

//test jakosci: sprawdzanie, ktory z nukeltoydow nie spelnia wymaganego progu jakosci
std::vector<std::vector<int>> qualityTest(std::vector<SEQ> instance, int qt) {
    std::vector<std::vector<int>> deleted(5);

    for (int i = 0 ; i < 5; i++ ) {
        for (int j = 0; j < instance[i].fasta.size(); j++) {
            if (instance[i].qual[j] < qt) {
                deleted[i].push_back(j);
            }
        }
    }

    return deleted;
}

//usuwanie z sekwencji nukletoydow o zbyt niskiej jakosci
std::vector<SEQ> clearInstace(std::vector <SEQ> instance, std::vector<std::vector<int>> toDelete){
    std::vector<SEQ> instance_checked(5);

    for(int i = 0; i < 5; i++){
        int t = 0;
        for(int j = 0; j < instance[i].fasta.size(); j++) {
            if(t > toDelete[i].size()) break;
            if(toDelete[i].empty()) goto copy;
            if (toDelete[i][t] == j){
                t++;
                continue;
            }
            copy:
                instance_checked[i].fasta.push_back(instance[i].fasta[j]);
                instance_checked[i].qual.push_back(instance[i].qual[j]);
                instance_checked[i].pos.push_back(instance[i].pos[j]);
        }
    }

    return instance_checked;
}

//tworzenie podciagow
std::vector<std::vector<std::vector<std::pair<char, int>>>> makeSubSeqs(int t, std::vector<SEQ> instance_checked){
    std::vector<std::vector<std::vector<std::pair<char, int>>>> instance_splited;

    for(int i  = 0; i < 5; i++){
        std::vector<std::vector<std::pair<char, int>>> subseqs;
        for(int j  = 0; j < instance_checked[i].fasta.size()-t+1; j++){
            std::vector<std::pair<char, int>> tmp;
            for(int k = 0; k < t+1; k++){
                std::pair<char, int> pair;
                pair = std::make_pair(instance_checked[i].fasta[j+k], instance_checked[i].pos[j+k]);
                tmp.push_back(pair);
            }
            subseqs.push_back(tmp);
            tmp.clear();
        }
        instance_splited.push_back(subseqs);
        subseqs.clear();
    }

    return instance_splited;
}

//łącznie w jeden rzad
std::vector<std::vector<std::pair<char, int>>> merge(std::vector<std::vector<std::vector<std::pair<char, int>>>> instance_splited){
    std::vector<std::vector<std::pair<char, int>>> merged;

    for(auto i: instance_splited){
        for(auto e: i){
            merged.push_back(e);
        }
    }
    return merged;
}

//generowanie grafu, gdzie kazdy wierzcholek to podciag, a krawdzie lacza takie same poociagi
int** generateGraph(std::vector<std::vector<std::pair<char, int>>> merged, int t, std::vector<std::vector<std::vector<std::pair<char, int>>>> instance_splited) {
    int n = merged.size();
    int subseq_size[6];
    subseq_size[0] = 0;
    for (int i = 0; i < 5; i++) {
        subseq_size[i + 1] = subseq_size[i] + instance_splited[i].size();
    }

    int **G = new int *[n];
    int u, v;
    for (u = 0; u < n; u++) {
        G[u] = new int[n];
    }
    for (u = 0; u < n; u++) {
        for (v = 0; v < n; v++) {
            G[u][v] = 0;
        }
    }

    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            int left = 0, right = 0;
            for(int x = 0; x < 5; x++){
                if(i > subseq_size[x] && i < subseq_size[x+1]){
                    left = i - subseq_size[x]+1;
                }
                if(j > subseq_size[x] && j < subseq_size[x+1]){
                    right = j - subseq_size[x]+1;
                }
            }
            if(left-right > t*10 || right-left > t*10) continue;
            int count = 0;
            for (int nuc = 0; nuc < t; nuc++) {
                if (merged[i][nuc].first == merged[j][nuc].first) count++;
            }
            if (count == t) {
                G[i][j] = 1;
                //G[j][i] = 1;
            }
        }
    }

    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            if (i == j) G[i][j] = 0;
        }
    }

    for (int i = 0; i < 5; i++) {
        for (int j = subseq_size[i]; j < subseq_size[i+1]; j++) {
            for (int m = subseq_size[i]; m < subseq_size[i+1]; m++) {
                G[j][m] = 0;
                //G[m][j] = 0;
            }
        }
    }

    return G;
}

//drukowanie grafu; testy
void printGraph(int** G, std::vector<std::vector<std::pair<char, int>>> merged){
    int n = merged.size();
    for (int k = 0; k < n; k++) {
        for (int j = 0; j < n; j++) {
            std::cout << G[k][j]<<" ";
        }
        std::cout << std::endl;
    }
}

//czy zawiera sie w wektorze
bool contains(std::vector<int> test, int x){
    for (auto e: test){
        if(x == e) return true;
    }
    return false;
}


//wyszukiwane struktury gwiazdy
std::vector<int> lookForStar(int** G, std::vector<std::vector<std::pair<char, int>>> merged){
    std::vector<int> cliques;

    for(int i = 0; i < merged.size(); i++){
        int s = 0;
        for(int j = 0; j < merged.size(); j++){
            if(G[i][j] == 1){
                s++;
            }
        }
        if(s > 3) cliques.push_back(i);
    }
    std::sort(cliques.begin(), cliques.end());
    return cliques;
}


//sprawdzanie motywu
std::vector<std::vector<int>> splitToTheme(std::vector<int> clique, std::vector<std::vector<std::pair<char, int>>> merged, std::vector<std::vector<std::vector<std::pair<char, int>>>> instance_splited, int t){
    int subseq_size[6];
    subseq_size[0] = 0;
    for (int i = 0; i < 5; i++) {
        subseq_size[i+1] = subseq_size[i] + instance_splited[i].size();
    }

    std::vector<std::vector<int>> theme;

    for(int i = 1; i < 6; i++){
        std::vector<int> tmp;
        for(auto e: clique){
            if(e > subseq_size[i-1]-1 && e < subseq_size[i]){
                tmp.push_back(e-subseq_size[i-1]);
            }
        }
        theme.push_back(tmp);
    }

    //usuwanie podwojnych???
    for(int i = 0; i < 5; i++){
        for(int j = 0; j < theme[i].size(); j++){
            for(int k = 0; k < theme[i].size(); k++){
                if(j == k) continue;
                int sum = 0;
                for(int m = 0; m < t; m++) {
                    if (instance_splited[i][theme[i][j]][m].first == instance_splited[i][theme[i][k]][m].first) sum++;
                }
                if(sum == t) theme[i].erase(theme[i].begin()+k);
            }
        }
    }

    for(int i = 0; i < 5; i++) {
        for (int o = 0; o < 5; o++) {
            if (i == o) continue;
            for (int j = 0; j < theme[i].size(); j++) {
                for (int k = 0; k < theme[o].size(); k++) {
                    int sum = 0;
                    for (int m = 0; m < t; m++) {
                        if (instance_splited[o][theme[o][k]][m].first == instance_splited[i][theme[i][j]][m].first) sum++;
                    }
                    if (sum == t) {
                        break;
                    }
                    if (k == theme[o].size() - 1) {
                        theme[i].erase(theme[i].begin() + j);
                        i = 0;
                        o = 0;
                    }
                }
            }
        }
    }

    return theme;
}

//drukowanie rezulatatu
void printTheme(std::vector<std::vector<int>> theme, std::vector<std::vector<std::vector<std::pair<char, int>>>> instance_splited, int t){

    std::cout<<"Znaleziono motyw w sekwencjach: ";

    for(int i = 0; i < 5; i++){
        std::cout<<std::endl<<"Sekwencja numer: "<< i+1<<std::endl;
        for(int cn = 0; cn < theme[i].size(); cn++){
            std::cout<<"podciag "<<cn+1<<" z "<< theme[i].size() <<": ";
            for(int j = 0; j < t; j++){
                std::cout<<instance_splited[i][theme[i][cn]][j].first;
            }
            std::cout<<" na pozycjach: ";
            for(int j = 0; j < t; j++){
                std::cout<<instance_splited[i][theme[i][cn]][j].second<<" ";
            }
            std::cout<<std::endl;
        }
    }
    std::cin.ignore();
}


//main
int main() {
    int t, qt;
    do {
        std::cout << "Podaj prog minimalnej wiarygodnosci (0-40): ";
        std::cin >> qt;
        if (qt < 0 || qt > 40) std::cout << "Wartosc musi byc z przedzialu 0-40!" << std::endl;
    } while (qt < 0 || qt > 40);
    do {
        std::cout << "Podaj dlugosc podciagow (4-9): ";
        std::cin >> t;
        if (t < 4 || t > 9) std::cout << "Wartosc musi byc z przedzialu 4-9!" << std::endl;
    } while (t < 4 || t > 9);

    std::vector<SEQ> instance(5);
    int instanceNum;
    std::cout << "Podaj numer instancji: ";
    std::cin >> instanceNum;
    std::cout << std::endl;
    std::cin.ignore();
    
    instance = getInstance(instanceNum);

    std::vector<std::vector<int>> toDelete;
    toDelete = qualityTest(instance, qt);

    std::vector<SEQ> instance_checked(5);
    instance_checked = clearInstace(instance, toDelete);

    for (auto e: instance_checked) {
        if (e.fasta.size() < t) {
            std::cout<< "Ustawiony zbyt wysoki prog oczekiwanej jakosci. Jedna z pieciu badanych sekwencji jest krotsza niz oczekiwana dlugosc podicagu.";
            std::cin.ignore();
            exit(3);
        }
    }

    std::vector<std::vector<std::vector<std::pair<char, int>>>> instance_splited;
    instance_splited = makeSubSeqs(t, instance_checked);

    std::vector<std::vector<std::pair<char, int>>> merged;
    merged = merge(instance_splited);

    int **G;
    G = generateGraph(merged, t, instance_splited);

    std::vector<int> clique;
    clique = lookForStar(G, merged);

        if(clique.size()  < 5){
        std::cout<<"Nie znaleziono motywu";
        std::cin.ignore();
        exit(4);
    }

    std::vector<std::vector<int>> theme;
    theme = splitToTheme(clique, merged, instance_splited, t);

    for(int i = 0; i < 5; i++){
        if(theme[i].empty()){
            std::cout<<"Nie znaleziono motywu";
            std::cin.ignore();
            
            exit(5);
        }
    }

    printTheme(theme, instance_splited, t);
	
    return 0;
}
