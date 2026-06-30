#ifndef HALIR_INTERNAL_H
#define HALIR_INTERNAL_H

#include <stddef.h>
#include <stdio.h>
#include <string.h>

static inline int
halir_copy_string_checked(char *dst, size_t dst_size, const char *src, const char *field_name)
{
  size_t src_len;

  if ((dst == NULL) || (src == NULL) || (dst_size == 0)) {
    fprintf(stderr, "Invalid string copy args for field: %s\n", field_name);
    return 0;
  }

  src_len = strlen(src);
  if (src_len >= dst_size) {
    fprintf(stderr, "Field '%s' exceeds max length (%zu)\n", field_name, dst_size - 1);
    return 0;
  }

  memcpy(dst, src, src_len + 1);
  return 1;
}

#endif