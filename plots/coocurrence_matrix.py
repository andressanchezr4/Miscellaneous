# Conteo absoluto de todas las veces que se da una pareja en una molecula
def co_occurrence(all_molecule_fg_list):
    co_dict = defaultdict(int)
    vocab = set()
    for fg_molecule_list in all_molecule_fg_list:
        window_size = len(fg_molecule_list)
        for i in range(len(fg_molecule_list)):
            token = fg_molecule_list[i]
            vocab.add(token)
            next_tokens = set(fg_molecule_list[i+1 : i+1+window_size]) #Lista con todos los tokens siguientes
            for t in next_tokens:
                key = tuple(sorted([t, token])) # Generamos la pair-key
                if key[1] == key[0]: # Filtramos las claves que sean dos grupos iguales
                    continue
                co_dict[key] += 1 # Añadimos una ocurrencia a la key
    
    # formulate the dictionary into dataframe
    vocab = sorted(vocab) # sort vocab
    co_matrix = pd.DataFrame(data=np.zeros((len(vocab), len(vocab)), dtype=np.int16),
                      index=vocab,
                      columns=vocab)
    for key, value in co_dict.items():
        co_matrix.at[key[0], key[1]] = value
        co_matrix.at[key[1], key[0]] = value
    return co_matrix, co_dict

c_matrix, c_dict = co_occurrence(fg_list)

# Conteo del numero de veces que se encuentra UNA pareja POR molecula
def co_occurrence_unique(all_molecule_fg_list):
    dictionary_list = []
    vocab = set()
    for fg_molecule_list in all_molecule_fg_list:
        window_size = len(fg_molecule_list)
        co_dict = defaultdict(int)
        for i in range(len(fg_molecule_list)):
            token = fg_molecule_list[i]
            vocab.add(token)
            next_tokens = set(fg_molecule_list[i+1 : i+1+window_size]) #Lista con todos los tokens siguientes
            for t in next_tokens:
                key = tuple(sorted([t, token])) # Generamos la pair-key
                if key[1] == key[0]: # Filtramos las claves que sean dos grupos iguales
                    continue
                co_dict[key] = 1 # Añadimos la primera ocurrencia a la key
        dictionary_list.append(co_dict)
    
    co_dict = {}
    for dictionary in dictionary_list:
        for key, values in dictionary.items():
            co_dict[key] = co_dict.get(key, 0) + 1
              
    # formulate the dictionary into dataframe
    vocab = sorted(vocab) # sort vocab
    co_matrix = pd.DataFrame(data=np.zeros((len(vocab), len(vocab)), dtype=np.int16),
                      index=vocab,
                      columns=vocab)
    for key, value in co_dict.items():
        co_matrix.at[key[0], key[1]] = value
        co_matrix.at[key[1], key[0]] = value
    return co_matrix, co_dict

c_unique_matrix, c_unique_dict = co_occurrence_unique(fg_list)

# Functional group pair unique ocurrence distribution
# Extraemos las X parejas que mas se repiten 
top100_fg_unique = Counter(c_unique_dict).most_common(100)
representacion_100 = sum([i[1] for i in top100_fg_unique])
representacion_total = sum(c_unique_dict.values())
representacion = representacion_100/representacion_total

top100_fg_unique_graph = top100_fg_unique[:30]
database = [str(fg[0][0] + ' - ' + fg[0][1]) for fg in top100_fg_unique_graph]
targets = [fg[1] for fg in top100_fg_unique_graph]
plt.barh(database, targets, color = ['darkorange'])
plt.title("Drugbank. Top 30 functional group pair ocurrences")
plt.ylabel("Functional groups", size = 12)
plt.xlabel('Frequency', size = 12)
plt.figure(figsize=(5,1))
plt.show()

# Redimensionamos la matriz de ocurrencia de forma que solo contenga información
# sobre las 40 parejas más frecuentes
vocab_0 = [i[0][0] for i in top100_fg_unique]
vocab_1 = [i[0][1] for i in top100_fg_unique]
vocab = set(vocab_1 + vocab_0)
co_matrix = pd.DataFrame(data=np.zeros((len(vocab), len(vocab)), dtype=np.int16),
                         index=vocab,
                         columns=vocab)

# Volvemos a hacer el diccionario con las 100 parejas más frecuentes y 
# construimos la matriz 
top100_fg_unique = {i[0]: i[1] for i in top100_fg_unique}
for key, value in top100_fg_unique.items():
        co_matrix.at[key[0], key[1]] = value
        co_matrix.at[key[1], key[0]] = value

