import csv
import itertools
import sys

PROBS = {

    # Unconditional probabilities for having gene
    "gene": {
        2: 0.01,
        1: 0.03,
        0: 0.96
    },

    "trait": {

        # Probability of trait given two copies of gene
        2: {
            True: 0.65,
            False: 0.35
        },

        # Probability of trait given one copy of gene
        1: {
            True: 0.56,
            False: 0.44
        },

        # Probability of trait given no gene
        0: {
            True: 0.01,
            False: 0.99
        }
    },

    # Mutation probability
    "mutation": 0.01
}


def main():

    # Check for proper usage
    if len(sys.argv) != 2:
        sys.exit("Usage: python heredity.py data.csv")
    people = load_data(sys.argv[1])

    # Keep track of gene and trait probabilities for each person
    probabilities = {
        person: {
            "gene": {
                2: 0,
                1: 0,
                0: 0
            },
            "trait": {
                True: 0,
                False: 0
            }
        }
        for person in people
    }

    # Loop over all sets of people who might have the trait
    names = set(people)
    for have_trait in powerset(names):

        # Check if current set of people violates known information
        fails_evidence = any(
            (people[person]["trait"] is not None and
             people[person]["trait"] != (person in have_trait))
            for person in names
        )
        if fails_evidence:
            continue

        # Loop over all sets of people who might have the gene
        for one_gene in powerset(names):
            for two_genes in powerset(names - one_gene):

                # Update probabilities with new joint probability
                p = joint_probability(people, one_gene, two_genes, have_trait)
                update(probabilities, one_gene, two_genes, have_trait, p)

    # Ensure probabilities sum to 1
    normalize(probabilities)

    # Print results
    for person in people:
        print(f"{person}:")
        for field in probabilities[person]:
            print(f"  {field.capitalize()}:")
            for value in probabilities[person][field]:
                p = probabilities[person][field][value]
                print(f"    {value}: {p:.4f}")


def load_data(filename):
    """
    Load gene and trait data from a file into a dictionary.
    File assumed to be a CSV containing fields name, mother, father, trait.
    mother, father must both be blank, or both be valid names in the CSV.
    trait should be 0 or 1 if trait is known, blank otherwise.
    """
    data = dict()
    with open(filename) as f:
        reader = csv.DictReader(f)
        for row in reader:
            name = row["name"]
            data[name] = {
                "name": name,
                "mother": row["mother"] or None,
                "father": row["father"] or None,
                "trait": (True if row["trait"] == "1" else
                          False if row["trait"] == "0" else None)
            }
    return data


def powerset(s):
    """
    Return a list of all possible subsets of set s.
    """
    s = list(s)
    return [
        set(s) for s in itertools.chain.from_iterable(
            itertools.combinations(s, r) for r in range(len(s) + 1)
        )
    ]


def joint_probability(people, one_gene, two_genes, have_trait):
    """
    Compute and return a joint probability.

    The probability returned should be the probability that
        * everyone in set `one_gene` has one copy of the gene, and
        * everyone in set `two_genes` has two copies of the gene, and
        * everyone not in `one_gene` or `two_gene` does not have the gene, and
        * everyone in set `have_trait` has the trait, and
        * everyone not in set` have_trait` does not have the trait.
    """
    

    #set for people with zero gene
    zero_gene = set()
    for person in people.keys():
        if person not in one_gene and person not in two_genes:
            zero_gene.add(person)

    #set for people with no trait
    no_trait = set()
    for person in people.keys():
        if person not in have_trait:
            no_trait.add(person)

    #calculate indivisual probability for those persons and their characteristics
    temp_prob = dict()
    for person in people.keys():
        temp_prob[person] = 0


    for person in people.keys():
        p = 1 
        mommy=people[person]['mother']
        daddy=people[person]['father']

        
        if person in one_gene:
            mommy=people[person]['mother']
            daddy=people[person]['father']
            if mommy==None and daddy == None:
                p = PROBS['gene'][1]
            else:
                #bad from mother and good from father
                if mommy in one_gene:
                    p1 = 0.5
                elif mommy in two_genes:
                    p1 = 1 - PROBS['mutation']
                elif mommy in zero_gene:
                    p1 = PROBS['mutation']

                if daddy in one_gene:
                    p1 *= 0.5
                elif daddy in two_genes:
                    p1 *= PROBS['mutation']
                elif daddy in zero_gene:
                    p1 *= 1 - PROBS['mutation']


                #good from mother and bad from father
                if mommy in one_gene:
                    p2 = 0.5
                elif mommy in two_genes:
                    p2 = PROBS['mutation']
                elif mommy in zero_gene:
                    p2 = 1 - PROBS['mutation']

                if daddy in one_gene:
                    p2*=0.5
                elif daddy in two_genes:
                    p2*= 1 - PROBS['mutation']
                elif daddy in zero_gene:
                    p2*= PROBS['mutation']

                p=p1+p2
                #p is now the probability that child has 1 gene given characters of parents
            p=p*PROBS['trait'][1][person in have_trait]
            #p is now the prob that states the trait also 
                        

        elif person in two_genes:

            mommy=people[person]['mother']
            daddy=people[person]['father']

            if mommy == None and daddy == None:
                p = PROBS['gene'][2]
            else:
                if mommy in one_gene:
                    p = 0.5
                elif mommy in two_genes:
                    p = 1 - PROBS['mutation']
                elif mommy in zero_gene:
                    p = PROBS['mutation']

                if daddy in one_gene:
                    p *= 0.5
                elif daddy in two_genes:
                    p *= 1 - PROBS['mutation']
                elif daddy in zero_gene:
                    p *= PROBS['mutation']

            p=p*PROBS['trait'][2][person in have_trait]


        elif person in zero_gene:
            if mommy== None and daddy == None:
                p = PROBS['gene'][0]
            else:
                if mommy in one_gene:
                    p = 0.5
                elif mommy in two_genes:
                    p=PROBS['mutation']
                elif mommy in zero_gene:
                    p=1 - PROBS['mutation']

                if daddy in one_gene:
                    p*= 0.5
                elif daddy in two_genes:
                    p*= PROBS['mutation']
                elif daddy in zero_gene:
                    p*= 1 - PROBS['mutation']

            p = p * PROBS['trait'][0][person in have_trait]

        temp_prob[person] = p


    cum_prob = 1
    for value in temp_prob.values():
        cum_prob *= value

    return cum_prob




        

        





def update(probabilities, one_gene, two_genes, have_trait, p):
    """
    Add to `probabilities` a new joint probability `p`.
    Each person should have their "gene" and "trait" distributions updated.
    Which value for each distribution is updated depends on whether
    the person is in `have_gene` and `have_trait`, respectively.
    """
    #raise NotImplementedError
    for person in probabilities:
        if person in one_gene:
            probabilities[person]['gene'][1] += p
        elif person in two_genes:
            probabilities[person]['gene'][2] += p
        else:
            probabilities[person]['gene'][0] += p

        if person in have_trait:
            probabilities[person]['trait'][True] += p
        else:
            probabilities[person]['trait'][False] += p


def normalize(probabilities):
    """
    Update `probabilities` such that each probability distribution
    is normalized (i.e., sums to 1, with relative proportions the same).
    """
    #raise NotImplementedError
    for person in probabilities:
        s1= sum(probabilities[person]["gene"].values())
        for key in probabilities[person]["gene"].keys():
            probabilities[person]["gene"][key] =probabilities[person]["gene"][key]/s1

        s2 = sum(probabilities[person]['trait'].values())
        for key in probabilities[person]['trait'].keys():
            probabilities[person]['trait'][key] = probabilities[person]['trait'][key]/s2

if __name__ == "__main__":
    main()
