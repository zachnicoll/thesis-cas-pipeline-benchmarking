class CasFamilyCount:
    actual_count: int
    true_positives: int
    false_positives: int

    def __init__(self) -> None:
        self.actual_count = 0
        self.true_positives = 0
        self.false_positives = 0

    def false_negatives(self) -> int:
        FNs = self.actual_count - (self.true_positives + self.false_positives)
        return FNs if FNs > 0 else 0  # Clamp count of FNs to 0
