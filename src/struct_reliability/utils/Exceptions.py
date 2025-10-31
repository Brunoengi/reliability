class InvalidKeyError(Exception):
    """Exception to invalid Keys."""
    def __init__(self, invalid_keys):
        super().__init__(f"Error: Invalid keys found: {invalid_keys}")
        self.invalid_keys = invalid_keys

