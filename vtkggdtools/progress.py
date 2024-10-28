class Progress:
    def __init__(self, update_progress):
        self.update_progress = update_progress
        self.progress = 0

    def increment(self, increment):
        """Increments the progress by a specified amount.

        Args:
            increment: amount to increment the progress by.
        """
        self.progress = self.progress + increment
        self.update_progress(self.progress)

    def set(self, value):
        """Sets the progress to a specific value (0 to 1).

        Args:
            value: specific value to set the progress to
        """
        self.progress = value
        self.update_progress(value)
