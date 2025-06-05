#include <iostream>
#include <algorithm>
#include <random>
#include <chrono>
#include <vector>
#include <fstream>

using namespace std;

struct Item
{
    int value, weight;
};

struct Indiv
{
    vector<int> genes;
    int score, weight;
};

const int POP_SIZE = 100;
const int GEN_COUNT = 1000;
const double MUT_RATE = 0.05;

random_device rd;
mt19937 gen(rd());

int evaluate(Indiv &ind, const vector<Item> &items, int maxWeight)
{
    int totalValue = 0, totalWeight = 0;
    for (size_t i = 0; i < ind.genes.size(); i++) {
        if (ind.genes[i]) {
            totalWeight += items[i].weight;
            totalValue += items[i].value;
        }
    }

    if (totalWeight > maxWeight) {
        return 0;
    }
    ind.weight = totalWeight;
    return totalValue;
}

Indiv getRandom(int size, const vector<Item> &items, int maxWeight)
{
    Indiv ind;
    ind.genes.resize(size, 0);
    ind.weight = 0;

    vector<int> indices(size);
    iota(indices.begin(), indices.end(), 0);
    shuffle(indices.begin(), indices.end(), gen);

    for (int i : indices) {
        if (ind.weight + items[i].weight <= maxWeight) {
            ind.genes[i] = 1;
            ind.weight += items[i].weight;
        }
    }
    ind.score = evaluate(ind, items, maxWeight);
    return ind;
}

void mutate(Indiv &ind, const vector<Item> &items, int maxWeight)
{
    uniform_real_distribution<double> mutationChance(0.0, 1.0);
    for (int i = 0; i < ind.genes.size(); i++) {
        if (mutationChance(gen) < MUT_RATE) {
            ind.genes[i] = 1 - ind.genes[i];
            if (evaluate(ind, items, maxWeight) == 0) {
                ind.genes[i] = 1 - ind.genes[i];
            }
        }
    }
    ind.score = evaluate(ind, items, maxWeight);
}

pair<Indiv, Indiv> crossOver(const Indiv &p1, const Indiv &p2, const vector<Item> &items, int maxWeight)
{
    int size = p1.genes.size();
    uniform_int_distribution<int> dist(0, size - 1);
    int cut = dist(gen);

    Indiv c1 = p1;
    Indiv c2 = p2;

    for (int i = cut; i < size; i++) {
        swap(c1.genes[i], c2.genes[i]);
    }
    c1.score = evaluate(c1, items, maxWeight);
    c2.score = evaluate(c2, items, maxWeight);

    if (c1.score == 0) {
        c1 = p1;
    }
    if (c2.score == 0) {
        c2 = p2;
    }

    return {c1, c2};
}

Indiv selectParent(const vector<Indiv> &pop)
{
    uniform_int_distribution<int> dist(0, pop.size() - 1);
    int a = dist(gen);
    int b = dist(gen);
    return (pop[a].score > pop[b].score) ? pop[a] : pop[b];
}

int SolveKS(vector<Item> &items, int maxWeight)
{
    int size = items.size();
    vector<Indiv> pop(POP_SIZE);
    for (Indiv &ind : pop) {
        ind = getRandom(size, items, maxWeight);
    }
    for (int i = 0; i < GEN_COUNT; i++) {
        vector<Indiv> newPop;
        while (newPop.size() < POP_SIZE) {
            Indiv p1 = selectParent(pop);
            Indiv p2 = selectParent(pop);
            auto [c1, c2] = crossOver(p1, p2, items, maxWeight);
            mutate(c1, items, maxWeight);
            mutate(c2, items, maxWeight);

            newPop.push_back(c1);
            if (newPop.size() < POP_SIZE) {
                newPop.push_back(c2);
            }
        }
        pop = newPop;
    }
    return max_element(pop.begin(), pop.end(), [](const Indiv &a, const Indiv &b) { return a.score < b.score;})->score;
}

void tester(string filename)
{
    auto totalStart = std::chrono::high_resolution_clock::now();
    cout << filename << endl;
    ifstream file(filename, ios::in);
    int n, w;
    file >> n >> w;
    vector<Item> items(n);
    for (auto &item : items) {
        file >> item.value >> item.weight;
    }
    file.close();
    auto algoStart =  std::chrono::high_resolution_clock::now();
    int result = SolveKS(items, w);
    auto algoEnd =  std::chrono::high_resolution_clock::now();
    auto totalEnd =  std::chrono::high_resolution_clock::now();
    cout << "result: " << result << endl;
    cout << "Algorithm time: " << std::chrono::duration_cast<std::chrono::microseconds>(algoEnd - algoStart).count() << " mks" << endl;
    cout << "Total time: " << std::chrono::duration_cast<std::chrono::microseconds>(totalEnd - totalStart).count() << " mks" << endl << endl;
    return;
}

int main() {
    tester("ks_4_0");
    tester("ks_30_0");
    tester("ks_100_0");
    tester("ks_300_0");
    return 0;
}
