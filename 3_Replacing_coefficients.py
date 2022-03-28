import pandas as pd
import re


def GetReactantsProductsStr(reaction):
    reactantproduct = reaction.split("<=>")
    reactantstr = reactantproduct[0]
    productstr = reactantproduct[1]
    return reactantstr, productstr


def GetReactantsProductsList(reactantstr, productstr):
    reactants = reactantstr.split(" + ")
    products = productstr.split(" + ")
    return reactants, products


def GetCoefficientsDict(metabolites):
    metdict = {}
    count = 0

    for metabolite in metabolites:
        coefficient = re.search(r"(?:^| )(\d[nmx]|[\dnmx]|\([nmx][+\-][\dmnx]\)|\([nmx]\))", str(metabolite))
        try:
            coefficient = coefficient.group(1)
        except AttributeError:
            coefficient = re.search(r"(\([nmx][+\-][\dmnx]\)|\([nmx]\))", str(metabolite))
            try:
                coefficient = coefficient.group(1)
            except AttributeError:
                coefficient = '1'

        compound = re.search(r"([CG]\d+)", str(metabolite))
        compound = compound.group(0) + "_" + str(count)

        metdict[compound] = coefficient
        count += 1
    return metdict


def GetMetsStr(metdict):
    metstr = ""

    for key, coefficient in metdict.items():
        compound = re.search(r'([CG]\d+)', str(key))
        compound = compound.group(0)
        if len(metstr) > 0:
            if coefficient == '1':
                metstr = metstr + " + " + compound
            else:
                metstr = metstr + " + " + coefficient + " " + compound
        else:
            if coefficient == '1':
                metstr = compound
            else:
                metstr = coefficient + " " + compound
    return metstr


def GetReactionStr(reactantdict, productdict):
    return GetMetsStr(metdict=reactantdict) + " <=> " + GetMetsStr(metdict=productdict)


# reactions_df = pd.read_csv("/home/masha/Research/Jupyter/KEGG/kegg_reactions_data_test.csv", sep=";")
# print("...Replacing coefficients for KEGG reactions...")
reactions_df = pd.read_csv("data/kegg_reactions.csv", sep=";")

for i in range(len(reactions_df)):
    reaction = reactions_df.Reaction[i]

    reactantstr, productstr = GetReactantsProductsStr(reaction=reaction)
    reactants, products = GetReactantsProductsList(reactantstr=reactantstr, productstr=productstr)
    reactantdict = GetCoefficientsDict(metabolites=reactants)
    # print("\n\nReactants: \n", reactantdict)
    productdict = GetCoefficientsDict(metabolites=products)
    # print("Products: \n", productdict, end="\n\n")

    # CASE 1
    # left: (n);
    # right: (n-x) AND (x) =>>
    # REPLACE (n) with 2, (x) with 1, (n-x) with 1 (R10976)

    if ("(n)" in reactantdict.values() or "n" in reactantdict.values()) & \
            ("(n-x)" in productdict.values() or "n-x" in productdict.values()) & \
            ("(x)" in productdict.values() or "x" in productdict.values()):

        for key, value in reactantdict.items():
            if value == 'n' or value == '(n)':
                reactantdict[key] = "2"

        for key, value in productdict.items():
            if value == 'x' or value == '(x)':
                productdict[key] = "1"
            elif value == '(n-x)' or value == 'n-x':
                productdict[key] = "1"

    # CASE 2
    # left: (m+n) OR (n+m);
    # right: ((m) AND (n)) OR (m AND n) OR m OR n =>>
    # replace (m+n) with 2, n with 1, (n) with 1, (m) with 1, m with 1

    # (m+n) on the left
    if ("(n+m)" in reactantdict.values() or "(m+n)" in reactantdict.values()) &\
            ((("(m)" in productdict.values() or "m" in productdict.values()) &
              ("(n)" in productdict.values() or "n" in productdict.values())) or
             ("(n)" in productdict.values() or "n" in productdict.values()) or
             ("(m)" in productdict.values() or "m" in productdict.values())):

        for key, value in reactantdict.items():
            if value == '(n+m)' or value == '(m+n)':
                reactantdict[key] = "2"

        for key, value in productdict.items():
            if value == 'm' or value == '(m)':
                productdict[key] = "1"
            elif value == 'n' or value == '(n)':
                productdict[key] = "1"

    # (m+n) on the right (mirrored version)
    elif ("(n+m)" in productdict.values() or "(m+n)" in productdict.values()) &\
            ((("(m)" in reactantdict.values() or "m" in reactantdict.values()) &
              ("(n)" in reactantdict.values() or "n" in reactantdict.values())) or
             ("(n)" in reactantdict.values() or "n" in reactantdict.values()) or
             ("(m)" in reactantdict.values() or "m" in reactantdict.values())):

        for key, value in productdict.items():
            if value == '(n+m)' or value == '(m+n)':
                reactantdict[key] = "2"

        for key, value in reactantdict.items():
            if value == 'm' or value == '(m)':
                productdict[key] = "1"
            elif value == 'n' or value == '(n)':
                productdict[key] = "1"

    # CASE 3
    # left: (m-1) AND (n+1);
    # right: (n OR(n)) AND (m OR (m)) =>>
    # replace (m-1) with 1, (m OR (m)) with 2, n+1 with 2, (n OR (n)) with 1 (R04193)

    if ("(n+1)" in reactantdict.values() or "n+1" in reactantdict.values()) & \
            ("(m-1)" in reactantdict.values() or "m-1" in reactantdict.values()) & \
            ("(n)" in productdict.values() or "n" in productdict.values()) & \
            ("(m)" in productdict.values() or "m" in productdict.values()):

        for key, value in reactantdict.items():
            if value == 'm-1' or value == '(m-1)':
                productdict[key] = "1"
            elif value == 'n+1' or value == '(n+1)':
                productdict[key] = "2"

        for key, value in productdict.items():
            if value == '(n)' or value == 'n':
                reactantdict[key] = "1"
            elif value == '(m)' or value == 'm':
                reactantdict[key] = "2"

    # mirrored version
    elif ("(n+1)" in productdict.values() or "n+1" in productdict.values()) & \
            ("(m-1)" in productdict.values() or "m-1" in productdict.values()) & \
            ("(n)" in reactantdict.values() or "n" in reactantdict.values()) & \
            ("(m)" in reactantdict.values() or "m" in reactantdict.values()):

        for key, value in productdict.items():
            if value == 'm-1' or value == '(m-1)':
                productdict[key] = "1"
            elif value == 'n+1' or value == '(n+1)':
                productdict[key] = "2"

        for key, value in reactantdict.items():
            if value == '(n)' or value == 'n':
                reactantdict[key] = "1"
            elif value == '(m)' or value == 'm':
                reactantdict[key] = "2"

    # CASE 4
    # left: (n-1) AND (m+1);
    # right: (n OR(n)) AND (m OR (m))=>>
    # replace (m+1) with 2, (m OR (m)) with 1, n-1 with 1, (n OR (n)) with 2 (R06184)

    if ("(n-1)" in reactantdict.values() or "n-1" in reactantdict.values()) & \
            ("(m+1)" in reactantdict.values() or "m+1" in reactantdict.values()) & \
            ("(n)" in productdict.values() or "n" in productdict.values()) & \
            ("(m)" in productdict.values() or "m" in productdict.values()):

        for key, value in reactantdict.items():
            if value == 'm+1' or value == '(m+1)':
                productdict[key] = "2"
            elif value == 'n-1' or value == '(n-1)':
                productdict[key] = "1"

        for key, value in productdict.items():
            if value == '(n)' or value == 'n':
                reactantdict[key] = "2"
            elif value == '(m)' or value == 'm':
                reactantdict[key] = "1"

    # mirrored version
    elif ("(n-1)" in productdict.values() or "n-1" in productdict.values()) & \
            ("(m+1)" in productdict.values() or "m+1" in productdict.values()) & \
            ("(n)" in reactantdict.values() or "n" in reactantdict.values()) & \
            ("(m)" in reactantdict.values() or "m" in reactantdict.values()):

        for key, value in productdict.items():
            if value == 'm+1' or value == '(m+1)':
                productdict[key] = "2"
            elif value == 'n-1' or value == '(n-1)':
                productdict[key] = "1"

        for key, value in reactantdict.items():
            if value == '(n)' or value == 'n':
                reactantdict[key] = "2"
            elif value == '(m)' or value == 'm':
                reactantdict[key] = "1"

    # CASE 5
    # left: n AND 2n AND 4n;
    # right: n AND 2n; =>>
    # replace 2n with 2, n with 1, 4n with 4

    if ("(n)" in reactantdict.values() or "n" in reactantdict.values()) & \
            ("(2n)" in reactantdict.values() or "2n" in reactantdict.values()) & \
            ("(4n)" in reactantdict.values() or "4n" in reactantdict.values()) & \
            ("(n)" in productdict.values() or "n" in productdict.values()) & \
            ("(2n)" in productdict.values() or "2n" in productdict.values()):

        for key, value in reactantdict.items():
            if value == 'n' or value == '(n)':
                reactantdict[key] = "1"
            elif value == '2n' or value == '(2n)':
                reactantdict[key] = "2"
            elif value == '4n' or value == '(4n)':
                reactantdict[key] = "4"

        for key, value in productdict.items():
            if value == 'n' or value == '(n)':
                productdict[key] = "1"
            elif value == '2n' or value == '(2n)':
                productdict[key] = "2"

    # CASE 6
    # left: n AND 2n;
    # right: n AND 2n AND (n+1); =>>
    # replace 2n with 2, n with 1, (n+1) with 2 (R05188, R05189)

    if ("(n)" in reactantdict.values() or "n" in reactantdict.values()) & \
            ("(2n)" in reactantdict.values() or "2n" in reactantdict.values()) & \
            ("(n)" in productdict.values() or "n" in productdict.values()) & \
            ("(2n)" in productdict.values() or "2n" in productdict.values()) & \
            ("(n+1)" in productdict.values() or "n+1" in productdict.values()):

        for key, value in reactantdict.items():
            if value == 'n' or value == '(n)':
                reactantdict[key] = "1"
            elif value == '2n' or value == '(2n)':
                reactantdict[key] = "2"

        for key, value in productdict.items():
            if value == 'n' or value == '(n)':
                productdict[key] = "1"
            elif value == '2n' or value == '(2n)':
                productdict[key] = "2"
            elif value == 'n+1' or value == '(n+1)':
                productdict[key] = "2"

    # CASE 7
    # left: n AND 2n;
    # right: n =>>
    # replace 2n with 2, n with 1 (R12328)

    if (("(n)" in reactantdict.values() or "n" in reactantdict.values()) &
        ("(2n)" in reactantdict.values() or "2n" in reactantdict.values()) &
        ("(n)" in productdict.values() or "n" in productdict.values())) or \
            (("(n)" in productdict.values() or "n" in productdict.values()) &
             ("(2n)" in productdict.values() or "2n" in productdict.values()) &
             ("(n)" in reactantdict.values() or "n" in reactantdict.values())):  # after "or" is a mirrored version

        for key, value in reactantdict.items():
            if value == 'n' or value == '(n)':
                reactantdict[key] = "1"
            elif value == '2n' or value == '(2n)':
                reactantdict[key] = "2"

        for key, value in productdict.items():
            if value == 'n' or value == '(n)':
                productdict[key] = "1"
            elif value == '2n' or value == '(2n)':
                productdict[key] = "2"

    # CASE 8
    # left: n;
    # right: 2n =>>
    # replace (2n) with 2, n with 1 (R05327)

    if (("(n)" in reactantdict.values() or "n" in reactantdict.values()) &
        ("(2n)" in productdict.values() or "2n" in productdict.values())) or \
            (("(n)" in productdict.values() or "n" in productdict.values()) &
             ("(2n)" in reactantdict.values() or "2n" in reactantdict.values())):  # after "or" is a mirrored version

        for key, value in reactantdict.items():
            if value == 'n' or value == '(n)':
                reactantdict[key] = "1"
            elif value == '2n' or value == '(2n)':
                reactantdict[key] = "2"

        for key, value in productdict.items():
            if value == 'n' or value == '(n)':
                productdict[key] = "1"
            elif value == '2n' or value == '(2n)':
                productdict[key] = "2"

    # CASE 9
    # left: ((n-2) OR n-2) AND (n OR (n));
    # right: (n-2) =>>
    # replace (n-2) with 1, n with 3, (n) with 3 (R10811)

    if (("(n-2)" in reactantdict.values() or "n-2" in reactantdict.values()) &
        ("(n)" in reactantdict.values() or "n" in reactantdict.values()) &
        ("(n-2)" in productdict.values() or "n-2" in productdict.values())) or \
            (("(n-2)" in productdict.values() or "n-2" in productdict.values()) &
             ("(n)" in productdict.values() or "n" in productdict.values()) &
             ("(n-2)" in reactantdict.values() or "n-2" in reactantdict.values())):  # after "or" is a mirrored version

        for key, value in reactantdict.items():
            if value == 'n' or value == '(n)':
                reactantdict[key] = "3"
            elif value == 'n-2' or value == '(n-2)':
                reactantdict[key] = "1"

        for key, value in productdict.items():
            if value == 'n' or value == '(n)':
                productdict[key] = "3"
            elif value == 'n-2' or value == '(n-2)':
                productdict[key] = "1"

    # CASE 10
    # left: ((n-1) OR n-1) AND (n OR (n));
    # right: n OR (n) =>>
    # replace (n-1) with 1, n with 2 (R10813)

    if (("(n-1)" in reactantdict.values() or "n-1" in reactantdict.values()) &
        ("(n)" in reactantdict.values() or "n" in reactantdict.values()) &
        ("(n)" in productdict.values() or "n" in productdict.values())) or \
            (("(n-1)" in productdict.values() or "n-1" in productdict.values()) &
             ("(n)" in productdict.values() or "n" in productdict.values()) &
             ("(n)" in reactantdict.values() or "n" in reactantdict.values())):  # after "or" is a mirrored version

        for key, value in reactantdict.items():
            if value == 'n' or value == '(n)':
                reactantdict[key] = "2"
            elif value == 'n-1' or value == '(n-1)':
                reactantdict[key] = "1"

        for key, value in productdict.items():
            if value == 'n' or value == '(n)':
                productdict[key] = "2"
            elif value == 'n-1' or value == '(n-1)':
                productdict[key] = "1"

    # CASE 11
    # left: (n-1) OR n-1
    # right: n OR (n) =>>
    # replace (n-1) with 1, n with 2 (R10813)

    if (("(n-1)" in reactantdict.values() or "n-1" in reactantdict.values()) &
        ("(n)" in productdict.values() or "n" in productdict.values())) or \
            (("(n-1)" in productdict.values() or "n-1" in productdict.values()) &
             ("(n)" in reactantdict.values() or "n" in reactantdict.values())):  # after "or" is a mirrored version

        for key, value in reactantdict.items():
            if value == 'n' or value == '(n)':
                reactantdict[key] = "2"
            elif value == 'n-1' or value == '(n-1)':
                reactantdict[key] = "1"

        for key, value in productdict.items():
            if value == 'n' or value == '(n)':
                productdict[key] = "2"
            elif value == 'n-1' or value == '(n-1)':
                productdict[key] = "1"

    # CASE 12
    # left: (n+1) OR n+1
    # right: n OR (n) =>>
    # replace (n+1) with 2, n with 1 (R10813)

    if (("(n+1)" in reactantdict.values() or "n+1" in reactantdict.values()) &
        ("(n)" in productdict.values() or "n" in productdict.values())) or \
            (("(n+1)" in productdict.values() or "n+1" in productdict.values()) &
             ("(n)" in reactantdict.values() or "n" in reactantdict.values())):  # after "or" is a mirrored version

        for key, value in reactantdict.items():
            if value == 'n' or value == '(n)':
                reactantdict[key] = "1"
            elif value == 'n+1' or value == '(n+1)':
                reactantdict[key] = "2"

        for key, value in productdict.items():
            if value == 'n' or value == '(n)':
                productdict[key] = "1"
            elif value == 'n+1' or value == '(n+1)':
                productdict[key] = "2"

    # CASE 13
    # left: (n+2) OR n+2
    # right: n OR (n) =>>
    # replace (n+2) with 2, n with 1 (R10813)

    if (("(n+2)" in reactantdict.values() or "n+2" in reactantdict.values()) &
        ("(n)" in productdict.values() or "n" in productdict.values())) or \
            (("(n+2)" in productdict.values() or "n+2" in productdict.values()) &
             ("(n)" in reactantdict.values() or "n" in reactantdict.values())):  # after "or" is a mirrored version

        for key, value in reactantdict.items():
            if value == 'n' or value == '(n)':
                reactantdict[key] = "1"
            elif value == 'n+2' or value == '(n+2)':
                reactantdict[key] = "3"

        for key, value in productdict.items():
            if value == 'n' or value == '(n)':
                productdict[key] = "1"
            elif value == 'n+2' or value == '(n+2)':
                productdict[key] = "3"

    # CASE 14
    # left: (n-1);
    # right: (n-1) =>>
    # replace (n-1) with 2

    if ("(n-1)" in reactantdict.values() or "n-1" in reactantdict.values()) & \
            ("(n-1)" in productdict.values() or "n-1" in productdict.values()):

        for key, value in reactantdict.items():
            if value == 'n-1' or value == '(n-1)':
                reactantdict[key] = "2"

        for key, value in productdict.items():
            if value == 'n-1' or value == '(n-1)':
                productdict[key] = "2"

    # CASE 15
    # left: (n-2);
    # right: (n-2) =>>
    # replace (n-2) with 1 (R02887)

    if ("(n-2)" in reactantdict.values() or "n-2" in reactantdict.values()) & \
            ("(n-2)" in productdict.values() or "n-2" in productdict.values()):

        for key, value in reactantdict.items():
            if value == 'n-2' or value == '(n-2)':
                reactantdict[key] = "3"

        for key, value in productdict.items():
            if value == 'n-2' or value == '(n-2)':
                productdict[key] = "3"

    # CASE 16
    # left: n OR (n);
    # right: n OR (n) =>>
    # replace n with 1, (n) with 1

    if ("(n)" in reactantdict.values() or "n" in reactantdict.values()) & \
            ("(n)" in productdict.values() or "n" in productdict.values()):

        for key, value in reactantdict.items():
            if value == 'n' or value == '(n)':
                reactantdict[key] = "1"

        for key, value in productdict.items():
            if value == 'n' or value == '(n)':
                productdict[key] = "1"

    # CASE 17
    # left: n OR (n) =>>
    # replace n with 1, (n) with 1

    if ("(n)" in reactantdict.values() or "n" in reactantdict.values()) or \
            ("(n)" in productdict.values() or "n" in productdict.values()):

        for key, value in reactantdict.items():
            if value == 'n' or value == '(n)':
                reactantdict[key] = "1"

        for key, value in productdict.items():
            if value == 'n' or value == '(n)':
                productdict[key] = "1"

    # print("Reactants:", reactantdict)
    # print("Products:", productdict)

    reaction = GetReactionStr(reactantdict, productdict)
    reactions_df.Reaction[i] = reaction


reactions_df.to_csv(r'data/kegg_reactions_redone.csv',
                    index=False)

reactions_df = reactions_df[["Reaction ID", "Reaction"]]
reactions_df.to_csv(r'data/kegg_reactions_redone.txt',
                    index=False, header=False, sep="\t")
