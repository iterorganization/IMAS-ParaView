class Progress:
    def __init__(self, update_progress=None, get_progress=None):
        self.get_progress = get_progress
        self.update_progress = update_progress

    def increment(self, increment):
        """Increments the progress by a specified amount.

        Args:
            increment: amount to increment the progress by.
        """
        if self.get_progress and self.update_progress:
            current_progress = self.get_progress()
            self.update_progress(current_progress + increment)

    def set(self, value):
        """Sets the progress to a specific value (0 to 1).

        Args:
            value: specific value to set the progress to
        """
        if self.update_progress:
            self.update_progress(value)
