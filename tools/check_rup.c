/* Forward RUP checker for comment-free generated DIMACS and addition-only proofs.
 * Every added clause is checked by unit propagation under its negation.
 * No SAT solver or clause-deletion implementation is part of this checker.
 */
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

typedef struct { int *literals; int length; } Clause;
static Clause *database;
static int count, capacity, variables;
static int8_t *assignment;

static void append(Clause clause) {
    if (count == capacity) {
        capacity = capacity ? 2 * capacity : 1024;
        database = realloc(database, (size_t)capacity * sizeof(Clause));
        if (!database) exit(2);
    }
    database[count++] = clause;
}

/* 1 = complete clause, 0 = clean EOF, -1 = malformed/truncated input. */
static int read_clause(FILE *input, Clause *clause) {
    int allocated = 16, length = 0, literal;
    int *values = malloc((size_t)allocated * sizeof(int));
    if (!values) return -1;
    for (;;) {
        int status = fscanf(input, "%d", &literal);
        if (status != 1) {
            free(values);
            return feof(input) && !length ? 0 : -1;
        }
        if (!literal) break;
        if (literal < -variables || literal > variables) {
            free(values);
            return -1;
        }
        if (length == allocated) {
            allocated *= 2;
            int *larger = realloc(values, (size_t)allocated * sizeof(int));
            if (!larger) { free(values); return -1; }
            values = larger;
        }
        values[length++] = literal;
    }
    *clause = (Clause){values, length};
    return 1;
}

static int set_literal(int literal) {
    int variable = literal > 0 ? literal : -literal;
    int value = literal > 0 ? 1 : -1;
    if (assignment[variable] == -value) return 0;
    assignment[variable] = (int8_t)value;
    return 1;
}

static int rup(Clause candidate) {
    memset(assignment, 0, (size_t)variables + 1);
    for (int i = 0; i < candidate.length; i++)
        if (!set_literal(-candidate.literals[i])) return 1;
    int changed;
    do {
        changed = 0;
        for (int i = 0; i < count; i++) {
            Clause clause = database[i];
            int remaining = 0, unit = 0, satisfied = 0;
            for (int j = 0; j < clause.length; j++) {
                int literal = clause.literals[j];
                int variable = literal > 0 ? literal : -literal;
                int value = assignment[variable];
                if (value == (literal > 0 ? 1 : -1)) { satisfied = 1; break; }
                if (!value) { remaining++; unit = literal; }
            }
            if (satisfied) continue;
            if (!remaining) return 1;
            if (remaining == 1) {
                if (!set_literal(unit)) return 1;
                changed = 1;
            }
        }
    } while (changed);
    return 0;
}

int main(int argc, char **argv) {
    if (argc != 2) return 2;
    FILE *input = fopen(argv[1], "r");
    if (!input) return 2;
    char prefix, format[4];
    int initial;
    if (fscanf(input, " %c %3s %d %d", &prefix, format, &variables, &initial) != 4 ||
        prefix != 'p' || strcmp(format, "cnf") || variables < 1 || initial < 0)
        return 2;
    assignment = calloc((size_t)variables + 1, sizeof(int8_t));
    if (!assignment) return 2;
    for (int i = 0; i < initial; i++) {
        Clause clause;
        if (read_clause(input, &clause) != 1) return 2;
        append(clause);
    }
    Clause trailing;
    if (read_clause(input, &trailing) != 0) return 2;
    fclose(input);
    int steps = 0, status;
    Clause clause;
    while ((status = read_clause(stdin, &clause)) == 1) {
        if (!rup(clause)) {
            fprintf(stderr, "RUP check failed at addition %d\n", steps + 1);
            return 1;
        }
        steps++;
        if (!clause.length) {
            printf("VERIFIED UNSAT additions %d variables %d initial_clauses %d\n",
                   steps, variables, initial);
            return 0;
        }
        append(clause);
    }
    fprintf(stderr, "Proof ended without a checked empty clause (status %d)\n", status);
    return 1;
}
