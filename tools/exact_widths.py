"""Integer widths for the proposed crossing construction, not a factoring engine."""


def seven_width(exponent: int) -> int:
    """Return ceil(log2(7**exponent)) without floating-point arithmetic."""
    if exponent < 0:
        raise ValueError("exponent must be nonnegative")
    return (7 ** exponent - 1).bit_length()


def piece_widths(label_width: int, index: int) -> tuple[int, int]:
    if label_width < 0:
        raise ValueError("label width must be nonnegative")
    blocks = (label_width + 2) // 3
    if not 0 <= index <= blocks:
        raise ValueError("piece index outside block range")
    return (
        min(2 * label_width, 4 * blocks + seven_width(index)),
        min(2 * label_width, 4 * blocks + seven_width(blocks - index)),
    )


def common_total_bound(grid_side: int, label_width: int) -> int:
    """Upper-bound the smaller ordered grid-network edge total."""
    if grid_side < 1 or label_width < 0:
        raise ValueError("require positive grid side and nonnegative label width")
    blocks = (label_width + 2) // 3
    width_sum = 8 * blocks + seven_width(blocks) + 1
    return grid_side ** 3 * label_width + grid_side ** 4 * width_sum // 4
