import re
from math import gcd
from functools import reduce
from sympy import Matrix


def parse_compound(compound):
    """
    Parses a chemical compound into a dictionary of elements and their counts.

    Args:
        compound (str): The chemical compound to parse.

    Returns:
        dict: A dictionary where keys are element symbols, and values are their counts.
    """
    elements = re.findall(r'([A-Z][a-z]*)(\d*)', compound)
    return {element: int(count) if count else 1 for element, count in elements}


def parse_reaction(reaction):
    """
    Parses the chemical reaction into reactants and products.

    Args:
        reaction (str): The reaction string to parse.

    Returns:
        tuple: Two lists, one for reactants and one for products.
    """
    try:
        left, right = reaction.split("->")
    except ValueError:
        return None, None
    
    reactants = [compound.strip() for compound in left.split("+")]
    products = [compound.strip() for compound in right.split("+")]
    
    return reactants, products


def get_elements(reactants, products):
    """
    Extracts all unique elements involved in the reaction from both reactants and products.

    Args:
        reactants (list): List of reactants (compounds).
        products (list): List of products (compounds).

    Returns:
        list: Sorted list of all unique elements in the reaction.
    """
    elements = set()
    for compound in reactants + products:
        elements.update(parse_compound(compound).keys())
    return sorted(elements)


def build_matrix(reactants, products, elements):
    """
    Constructs the matrix of element balances for the reaction.

    Args:
        reactants (list): List of reactant compounds.
        products (list): List of product compounds.
        elements (list): List of unique elements involved in the reaction.

    Returns:
        list: Matrix representing the coefficients of the elements in the compounds.
    """
    matrix = []
    for element in elements:
        row = []
        for compound in reactants:
            row.append(parse_compound(compound).get(element, 0))
        for compound in products:
            row.append(-parse_compound(compound).get(element, 0))
        matrix.append(row)
    return matrix


def simplify_coefficients(coefficients):
    """
    Simplifies the list of coefficients to their smallest whole numbers.

    Args:
        coefficients (list): List of coefficients to simplify.

    Returns:
        list: Simplified coefficients as integers.
    """
    def lcm(a, b):
        """Calculate the least common multiple of two numbers."""
        return abs(a * b) // gcd(a, b)
    
    denominators = [abs(f.denominator) for f in coefficients]
    scale_factor = reduce(lcm, denominators, 1)
    
    scaled_coeffs = [int(c * scale_factor) for c in coefficients]
    gcd_coeffs = reduce(gcd, scaled_coeffs)
    return [c // gcd_coeffs for c in scaled_coeffs]


def balance_reaction(reaction):
    """
    Balances a chemical reaction using a systematic method.

    Args:
        reaction (str): The chemical reaction to balance in the form "reactants -> products".

    Returns:
        str: The balanced chemical reaction.
    """
    reactants, products = parse_reaction(reaction)
    if reactants is None or products is None:
        return "Invalid reaction format. Use '->' to separate reactants and products."
    
    elements = get_elements(reactants, products)
    matrix = build_matrix(reactants, products, elements)
    
    m = Matrix(matrix)
    null_space = m.nullspace()
    if not null_space:
        return "This reaction cannot be balanced."
    
    coefficients = null_space[0]
    coefficients = simplify_coefficients(coefficients)
    
    # Separate reactant and product coefficients
    num_reactants = len(reactants)
    reactant_coeffs = coefficients[:num_reactants]
    product_coeffs = coefficients[num_reactants:]
    
    # Build the balanced reaction string
    balanced_reactants = " + ".join(f"{reactant_coeffs[i]} {reactants[i]}" for i in range(num_reactants))
    balanced_products = " + ".join(f"{product_coeffs[i]} {products[i]}" for i in range(len(products)))
    
    return f"{balanced_reactants} -> {balanced_products}"


if __name__ == "__main__":
    """
    The main function for running the reaction balancer in an interactive mode.
    It continuously asks for reactions to balance until the user decides to stop.
    """
    print("Chemical Reaction Balancer")
    print("Enter a chemical reaction to balance (e.g., H2 + O2 -> H2O):")

    while True:
        reaction = input("Reaction: ").strip()
        if reaction.lower() in ["exit", "quit"]:
            print("Program stopped.")
            break
        result = balance_reaction(reaction)
        print(f"Balanced Reaction: {result}")
