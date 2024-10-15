import sys


def parse_grammar_from_file(filename):
    """Parses a grammar from a file and returns it as a dictionary.

    Args:
        filename (str): The name of the file containing the grammar.

    Returns:
        dict: A dictionary representing the grammar.
    """
    grammar = {}
    current_non_terminal = None

    with open(filename, "r") as file:
        lines = file.readlines()

    for line in lines:
        line = line.strip()
        if not line:
            continue

        # Si es un nuevo no terminal
        if ":" in line:
            non_terminal, productions = line.split(":")
            non_terminal = non_terminal.strip()

            # Si el no terminal ya existe, no lo sobreescribimos, solo
            # agregamos producciones nuevas
            if non_terminal not in grammar:
                grammar[non_terminal] = []

            # Si hay producciones después de ':', agregarlas
            # (cuando no hay '|')
            if productions.strip():
                grammar[non_terminal].append(productions.strip().split())
            current_non_terminal = non_terminal

        # Si hay producciones separadas por '|', se añaden debajo del
        # no terminal
        elif "|" in line and current_non_terminal:
            for prod in line.split("|"):
                grammar[current_non_terminal].append(prod.strip().split())

    # Eliminar listas vacías en las producciones
    for non_terminal in grammar:
        grammar[non_terminal] = [prod for prod in grammar[non_terminal]
                                 if prod]

    return grammar


def extract_terminals_and_non_terminals(grammar):
    """Extracts terminals and non-terminals from the grammar.

    Args:
        grammar (dict): The grammar from which to extract terminals and
        non-terminals.

    Returns:
        tuple: A tuple containing a list of terminals and a list of
        non-terminals.
    """
    terminals = set()
    non_terminals = set(grammar.keys())

    for productions in grammar.values():
        for production in productions:
            for symbol in production:
                # Agregar a terminales si no es un no terminal
                if symbol not in non_terminals and symbol != "ε":
                    # Considerar letras minúsculas o símbolos no alfabéticos
                    # como terminales
                    if symbol.islower() or not symbol[0].isalpha():
                        terminals.add(symbol)

    return list(terminals), list(non_terminals)


def removeLeftRecursion(rulesgrammar):
    """Removes left recursion from the grammar rules.

    Args:
        rulesgrammar (dict): The grammar rules from which to remove left
        recursion.

    Returns:
        dict: The grammar rules with left recursion removed.
    """
    store = {}
    for lhs in rulesgrammar:
        alphaRules = []
        betaRules = []
        allrhs = rulesgrammar[lhs]
        for subrhs in allrhs:
            if subrhs[0] == lhs:
                alphaRules.append(subrhs[1:])
            else:
                betaRules.append(subrhs)
        if len(alphaRules) != 0:
            lhs_ = lhs + "'"
            while (lhs_ in rulesgrammar.keys()) or (lhs_ in store.keys()):
                lhs_ += "'"
            for b in range(0, len(betaRules)):
                betaRules[b].append(lhs_)
            rulesgrammar[lhs] = betaRules
            for a in range(0, len(alphaRules)):
                alphaRules[a].append(lhs_)
            alphaRules.append(["#"])
            store[lhs_] = alphaRules
    for left in store:
        rulesgrammar[left] = store[left]
    return rulesgrammar


def LeftFactoring(rulesgrammar):
    """Applies left factoring to the grammar rules.

    Args:
        rulesgrammar (dict): The grammar rules to which left factoring will be
        applied.

    Returns:
        dict: The grammar rules after left factoring.
    """
    newDict = {}
    for lhs in rulesgrammar:
        allrhs = rulesgrammar[lhs]
        temp = dict()
        for subrhs in allrhs:
            if subrhs[0] not in list(temp.keys()):
                temp[subrhs[0]] = [subrhs]
            else:
                temp[subrhs[0]].append(subrhs)
        new_rule = []
        tempo_dict = {}
        for term_key in temp:
            allStartingWithTermKey = temp[term_key]
            if len(allStartingWithTermKey) > 1:
                lhs_ = lhs + "'"
                while (lhs_ in rulesgrammar.keys()) or \
                      (lhs_ in tempo_dict.keys()):
                    lhs_ += "'"
                new_rule.append([term_key, lhs_])
                ex_rules = []
                for g in temp[term_key]:
                    ex_rules.append(g[1:])
                tempo_dict[lhs_] = ex_rules
            else:
                new_rule.append(allStartingWithTermKey[0])
        newDict[lhs] = new_rule
        for key in tempo_dict:
            newDict[key] = tempo_dict[key]
    return newDict


def save_grammar_to_file(grammar, output_filename):
    """Saves the grammar to a specified output file.

    Args:
        grammar (dict): The grammar to save.
        output_filename (str): The name of the output file.
    """
    with open(output_filename, "w") as file:
        for non_terminal, productions in grammar.items():
            productions_str = []
            for production in productions:
                if not production:  # Manejar el caso de ε
                    productions_str.append("ε")
                else:
                    productions_str.append(" ".join(filter(None, production)))

            # Formatear la línea de salida
            if len(productions_str) == 1:
                file.write(f"{non_terminal}: {productions_str[0]}\n")
            else:
                file.write(f"{non_terminal}:\n")
                for prod in productions_str:
                    file.write(f"    | {prod}\n")


if __name__ == "__main__":
    # Inicializar variables globales
    grammar = {}
    firsts = {}
    follows = {}
    term = []
    nonterm = []
    start_symbol = ""
    output_file = "gramatica_transformada.txt"

    # Comprobar si se pasa el archivo de gramática como argumento
    if len(sys.argv) < 2:
        print("Por favor, proporciona un archivo de gramática.")
        sys.exit(1)

    # Leer la gramática desde el archivo
    grammar = parse_grammar_from_file(sys.argv[1])
    start_symbol = list(grammar.keys())[0]

    print("\nRules: \n")
    for y in grammar:
        print(f"{y}->{grammar[y]}")

    print("\nAfter elimination of left recursion:\n")
    grammar = removeLeftRecursion(grammar)
    for y in grammar:
        print(f"{y}->{grammar[y]}")

    print("\nAfter left factoring:\n")
    grammar = LeftFactoring(grammar)
    for y in grammar:
        print(f"{y}->{grammar[y]}")

    term, nonterm = extract_terminals_and_non_terminals(grammar)
    print("\nTerminales: ", term)
    print("No terminales: ", nonterm)

    # Eliminar recursión izquierda
    grammar = removeLeftRecursion(grammar)
    # Aplicar factorización izquierda
    grammar = LeftFactoring(grammar)

    # Guardar gramática después de transformaciones
    save_grammar_to_file(grammar, output_file)
