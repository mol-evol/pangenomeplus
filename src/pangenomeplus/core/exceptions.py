"""Custom exceptions for the pangenome pipeline."""

import time
from typing import Optional, Dict, Any, List
from pathlib import Path


class PipelineError(Exception):
    """Base exception for pipeline errors."""

    def __init__(self, message: str, stage: Optional[str] = None) -> None:
        self.stage = stage
        self.timestamp = time.time()
        super().__init__(message)


class ToolExecutionError(PipelineError):
    """External tool execution failed."""

    def __init__(self, tool: str, command: str, returncode: int, stderr: str,
                 stage: Optional[str] = None) -> None:
        self.tool = tool
        self.command = command
        self.returncode = returncode
        self.stderr = stderr
        super().__init__(f"{tool} failed with return code {returncode}", stage)

    def get_error_details(self) -> Dict[str, Any]:
        """Get structured error details for logging/recovery."""
        return {
            "tool": self.tool,
            "command": self.command,
            "returncode": self.returncode,
            "stderr": self.stderr,
            "timestamp": self.timestamp,
            "stage": self.stage
        }


class ValidationError(PipelineError):
    """Data validation failed."""

    def __init__(self, message: str, errors: Optional[List[str]] = None,
                 stage: Optional[str] = None) -> None:
        self.errors = errors or []
        super().__init__(message, stage)


class ConfigurationError(PipelineError):
    """Configuration error."""

    def __init__(self, message: str, config_path: Optional[Path] = None,
                 stage: Optional[str] = None) -> None:
        self.config_path = config_path
        super().__init__(message, stage)


class ResourceError(PipelineError):
    """Insufficient system resources."""

    def __init__(self, message: str, resource_type: str, required: Optional[str] = None,
                 available: Optional[str] = None, stage: Optional[str] = None) -> None:
        self.resource_type = resource_type
        self.required = required
        self.available = available
        super().__init__(message, stage)


class DependencyError(PipelineError):
    """Missing or incompatible dependency."""

    def __init__(self, message: str, dependency: str, version_required: Optional[str] = None,
                 version_found: Optional[str] = None, stage: Optional[str] = None) -> None:
        self.dependency = dependency
        self.version_required = version_required
        self.version_found = version_found
        super().__init__(message, stage)