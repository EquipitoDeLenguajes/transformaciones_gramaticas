import sys
import copy


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


def first(rule):
    """Calculates the FIRST set for a given rule.

    Args:
        rule (list): The rule for which to calculate the FIRST set.

    Returns:
        list or str: The FIRST set of the rule.
    """
    global firsts, grammar, term
    if len(rule) != 0 and (rule is not None):
        if rule[0] in term:
            return rule[0]
        elif rule[0] == "#":
            return "#"

    if len(rule) != 0:
        if rule[0] in list(grammar.keys()):
            fres = []
            rhs_rules = grammar[rule[0]]
            for itr in rhs_rules:
                indivRes = first(itr)
                if type(indivRes) is list:
                    for i in indivRes:
                        fres.append(i)
                else:
                    fres.append(indivRes)

            if "#" not in fres:
                return fres
            else:
                newList = []
                fres.remove("#")
                if len(rule) > 1:
                    ansNew = first(rule[1:])
                    if ansNew is not None:
                        if type(ansNew) is list:
                            newList = fres + ansNew
                        else:
                            newList = fres + [ansNew]
                    else:
                        newList = fres
                    return newList
                fres.append("#")
                return fres


def follow(nt):
    """Calculates the FOLLOW set for a given non-terminal.

    Args:
        nt (str): The non-terminal for which to calculate the FOLLOW set.

    Returns:
        list: The FOLLOW set of the non-terminal.
    """
    global start_symbol, grammar
    solset = set()
    if nt == start_symbol:
        solset.add("$")

    for curNT in grammar:
        rhs = grammar[curNT]
        for subrule in rhs:
            if nt in subrule:
                while nt in subrule:
                    index_nt = subrule.index(nt)
                    subrule = subrule[index_nt + 1:]
                    if len(subrule) != 0:
                        res = first(subrule)
                        if "#" in res:
                            newList = []
                            res.remove("#")
                            ansNew = follow(curNT)
                            if ansNew is not None:
                                if type(ansNew) is list:
                                    newList = res + ansNew
                                else:
                                    newList = res + [ansNew]
                            else:
                                newList = res
                            res = newList
                    else:
                        if nt != curNT:
                            res = follow(curNT)

                    if res is not None:
                        if type(res) is list:
                            for g in res:
                                solset.add(g)
                        else:
                            solset.add(res)
    return list(solset)


def computeAllFirsts():
    """Computes the FIRST sets for all non-terminals in the grammar."""
    global grammar, term, nonterm

    for y in list(grammar.keys()):
        t = set()
        for sub in grammar.get(y):
            res = first(sub)
            if res is not None:
                if type(res) is list:
                    for u in res:
                        t.add(u)
                else:
                    t.add(res)

        firsts[y] = t

    print("\nCalculated firsts: ")
    key_list = list(firsts.keys())
    index = 0
    for gg in firsts:
        print(f"first({key_list[index]}) " f"=> {firsts.get(gg)}")
        index += 1


def computeAllFollows():
    """Computes the FOLLOW sets for all non-terminals in the grammar."""
    global grammar, follows
    for NT in grammar:
        solset = set()
        sol = follow(NT)
        if sol is not None:
            for g in sol:
                solset.add(g)
        follows[NT] = solset

    print("\nCalculated follows: ")
    key_list = list(follows.keys())
    index = 0
    for gg in follows:
        print(f"follow({key_list[index]})" f" => {follows[gg]}")
        index += 1


def createParseTable():
    """Creates the parsing table for the grammar and checks if it is LL(1).

    Returns:
        bool: True if the grammar is LL(1), False otherwise.
    """
    global grammar, firsts, follows, term
    print("\nFirsts and Follow Result table\n")

    # find space size
    mx_len_first = 0
    mx_len_fol = 0
    for u in grammar:
        k1 = len(str(firsts[u]))
        k2 = len(str(follows[u]))
        if k1 > mx_len_first:
            mx_len_first = k1
        if k2 > mx_len_fol:
            mx_len_fol = k2

    print(
        f"{{:<{10}}} "
        f"{{:<{mx_len_first + 5}}} "
        f"{{:<{mx_len_fol + 5}}}".format("Non-T", "FIRST", "FOLLOW")
    )
    for u in grammar:
        print(
            f"{{:<{10}}} "
            f"{{:<{mx_len_first + 5}}} "
            f"{{:<{mx_len_fol + 5}}}".format(u, str(firsts[u]),
                                             str(follows[u]))
        )

    # create matrix of row(NT) x [col(T) + 1($)]
    # create list of non-terminals
    ntlist = list(grammar.keys())
    terminals = copy.deepcopy(term)
    terminals.append("$")

    # create the initial empty state of ,matrix
    mat = []
    for x in grammar:
        row = []
        for y in terminals:
            row.append("")
        # of $ append one more col
        mat.append(row)

    # Classifying grammar as LL(1) or not LL(1)
    grammar_is_LL = True

    # rules implementation
    for lhs in grammar:
        rhs = grammar[lhs]
        for y in rhs:
            res = first(y)
            # epsilon is present,
            # - take union with follow
            if "#" in res:
                if type(res) is str:
                    firstFollow = []
                    fol_op = follows[lhs]
                    if fol_op is str:
                        firstFollow.append(fol_op)
                    else:
                        for u in fol_op:
                            firstFollow.append(u)
                    res = firstFollow
                else:
                    res.remove("#")
                    res = list(res) + list(follows[lhs])
            # add rules to table
            ttemp = []
            if type(res) is str:
                ttemp.append(res)
                res = copy.deepcopy(ttemp)
            for c in res:
                xnt = ntlist.index(lhs)
                yt = terminals.index(c)
                if mat[xnt][yt] == "":
                    mat[xnt][yt] = mat[xnt][yt] + f"{lhs}->{' '.join(y)}"
                else:
                    # if rule already present
                    if f"{lhs}->{y}" in mat[xnt][yt]:
                        continue
                    else:
                        grammar_is_LL = False
                        mat[xnt][yt] = mat[xnt][yt] + f",{lhs}->{' '.join(y)}"

    # final state of parse table
    print("\nGenerated parsing table:\n")
    frmt = "{:>12}" * len(terminals)
    print(frmt.format(*terminals))

    j = 0
    for y in mat:
        frmt1 = "{:>12}" * len(y)
        print(f"{ntlist[j]} {frmt1.format(*y)}")
        j += 1

    return (mat, grammar_is_LL, terminals)


if __name__ == "__main__":
    # Inicializar variables globales
    grammar = {}
    firsts = {}
    follows = {}
    term = []
    nonterm = []
    start_symbol = ""
    output_file = "output.txt"

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

    # Calcular conjuntos FIRST y FOLLOW
    computeAllFirsts()
    computeAllFollows()

    # Crear tabla de análisis
    is_ll1 = createParseTable()
    if is_ll1:
        print("La gramática es LL(1).")
    else:
        print("La gramática NO es LL(1).")
