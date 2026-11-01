"""Eliminate fresh odd relation pivots while preserving exact group identities.

The caller supplies genuine, rationally independent integer relation
rows. The engine retains a common odd scale D and integer words V_p
with D*g_p = V_p in the full group. It never cancels an odd common
factor from an expression. Its presented two-primary torsion is the
same as that of the uncompressed relation presentation.
"""


def _add(destination, source, multiplier):
    for label, value in source.items():
        destination[label] = destination.get(label, 0) + multiplier * value
        if destination[label] == 0:
            del destination[label]


class OddRelationPresentation:
    """Stream independent raw relations indexed by distinct integer labels.

    ``initial_labels`` name known group elements. A label first appearing
    in a later row is fresh. The largest fresh label with odd coefficient
    is eliminated; all other fresh labels become active generators.

    Neither group equality nor rational row independence is tested here.
    A caller must establish both. In the leading-norm application, the
    original distinct weighted-dominant pivot columns establish rank.
    """

    def __init__(self, initial_labels=()):
        labels = list(initial_labels)
        if (any(not isinstance(label, int) for label in labels)
                or len(set(labels)) != len(labels)):
            raise ValueError('Initial labels must be distinct integers.')
        self.scale = 1
        self.active_labels = labels
        self.expressions = {label: {label: 1} for label in labels}
        self.rows = []
        self.steps = []

    def add_relation(self, relation):
        """Add one genuine independent row; return its elimination record."""
        if (not isinstance(relation, dict)
                or any(not isinstance(label, int) or not isinstance(value, int)
                       for label, value in relation.items())):
            raise ValueError('A raw relation must map integer labels to integer entries.')
        relation = {label: value for label, value in relation.items() if value}
        if not relation:
            raise ValueError('An independent raw relation cannot be zero.')
        fresh = sorted(set(relation) - set(self.expressions))
        odd = [label for label in fresh if relation[label] % 2]
        eliminated = max(odd) if odd else None
        old_scale = self.scale
        for label in fresh:
            if label != eliminated:
                self.active_labels.append(label)
                self.expressions[label] = {label: old_scale}
        row = {}
        for label, coefficient in relation.items():
            if label != eliminated:
                _add(row, self.expressions[label], coefficient)
        pivot = None
        if eliminated is not None:
            pivot = relation[eliminated]
            multiplier = abs(pivot)
            sign = 1 if pivot > 0 else -1
            # W + c*D*g_q = 0. Multiply every older expression by |c|,
            # but set the new one to -sign(c)*W without cancelling D.
            for label in self.expressions:
                self.expressions[label] = {
                    active: multiplier * value
                    for active, value in self.expressions[label].items()}
            self.expressions[eliminated] = {
                active: -sign * value for active, value in row.items()}
            self.scale *= multiplier
        else:
            if not row:
                raise ValueError('The supplied raw independence promise failed.')
            self.rows.append(row)
            if len(self.rows) > len(self.active_labels):
                raise ValueError('There are more retained rows than active generators.')
        assert self.scale % 2 == 1
        record = {
            'raw_relation': relation, 'fresh_labels': fresh,
            'eliminated_label': eliminated, 'odd_pivot': pivot,
            'old_scale': old_scale, 'new_scale': self.scale,
            'substituted_relation': row, 'stored_relation': eliminated is None,
            'active_generators': len(self.active_labels),
            'retained_rows': len(self.rows),
        }
        self.steps.append(record)
        return record

    def retained_matrix(self):
        return [[row.get(label, 0) for label in self.active_labels] for row in self.rows]

    def snapshot(self):
        return {'odd_scale': self.scale, 'active_labels': list(self.active_labels),
                'expressions': {label: dict(word) for label, word in self.expressions.items()},
                'relation_rows': [dict(row) for row in self.rows],
                'steps': list(self.steps)}
